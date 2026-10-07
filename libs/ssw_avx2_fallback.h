/* Exact affine SW safety path, used only on SIMD overflow/invalid legacy
 * traceback. Scores use 32 bits; traceback uses one byte per DP cell.
 * The public SSW uint16_t score ABI cannot represent scores above 65535. */
#include <errno.h>
#include <limits.h>

static s_align* ssw_avx2_scalar(const s_profile* p, const int8_t* ref, int nr,
        uint8_t go, uint8_t ge, uint8_t flag, uint16_t filters, int filterd, int mask) {
    const int nq=p->readLen;
    const size_t stride=(size_t)nq+1;
    int *h=NULL,*e=NULL,*columns=NULL,*mstate=NULL,*fstate=NULL;
    uint8_t *dir=NULL;
    s_align* a=(s_align*)calloc(1,sizeof(*a));
    if(!a) return NULL;
    if(flag & SSW_REPORT_PREFIX_SCORES) {
        a->prefix_scores=(uint16_t*)calloc(nq,sizeof(uint16_t));
        if(!a->prefix_scores) {free(a);errno=ENOMEM;return NULL;}
    }
    flag &= ~SSW_REPORT_PREFIX_SCORES;
    a->ref_begin1=a->read_begin1=-1;
    h=(int*)calloc(stride,sizeof(int));
    e=(int*)calloc(stride,sizeof(int));
    /* When open < extend, reopening the SAME gap from H is not a valid
     * affine transition. Keep M/F separately for this unusual scoring. */
    if(go<ge) {mstate=(int*)calloc(stride,sizeof(int));fstate=(int*)calloc(stride,sizeof(int));}
    if(mask>=15) columns=(int*)calloc(nr,sizeof(int));
    if(flag) {
        if((size_t)nr>SIZE_MAX/(size_t)nq) { errno=ENOMEM; goto error; }
        dir=(uint8_t*)malloc((size_t)nr*nq);
    }
    if(!h||!e||(mask>=15&&!columns)||(flag&&!dir)||(go<ge&&(!mstate||!fstate))) { errno=ENOMEM; goto error; }
    for(int i=0;i<nr;++i) {
        int diagonal=0,f=0,column=0,ml=0,el=0;
        for(int j=1;j<=nq;++j) {
            int old=h[j],eo=old-go,ee=e[j]-ge,fo=h[j-1]-go,fe=f-ge;
            uint8_t d=0,ep=0,fp=0;
            if(mstate) {
                eo=(mstate[j]>=fstate[j]?mstate[j]:fstate[j])-go;
                ep=mstate[j]>=fstate[j]?1:3;
                fo=(ml>=el?ml:el)-go;fp=ml>=el?1:2;
            }
            if(eo<=ee)ep=2;
            if(fo<=fe)fp=3;
            e[j]=eo>ee?eo:ee; f=fo>fe?fo:fe;
            int v=diagonal+p->mat[ref[i]*p->n+p->read[j-1]];
            if(v>0) d=1; else v=0;
            if(a->prefix_scores && v>a->prefix_scores[j-1]) a->prefix_scores[j-1]=v;
            if(mstate) {mstate[j]=ml=v;fstate[j]=f;el=e[j];}
            if(e[j]>v) {v=e[j];d=2;}
            if(f>v) {v=f;d=3;}
            if(dir) dir[(size_t)i*nq+j-1]=d|(ep<<2)|(fp<<4);
            h[j]=v;diagonal=old;
            if(v>65535) {errno=ERANGE;goto error;}
            if(v>column) column=v;
            if(v>a->score1) {a->score1=v;a->ref_end1=i;a->read_end1=j-1;}
        }
        if(columns) columns[i]=column;
    }
    if(!a->score1) goto done;
    a->ref_end2=mask>=15?0:-1;
    if(columns) for(int i=0;i<nr;++i) {
        if((i<a->ref_end1-mask||i>=a->ref_end1+mask)&&columns[i]>a->score2) {
            a->score2=columns[i];a->ref_end2=i;
        }
    }
    if(!flag||(flag==2&&a->score1<filters)) goto done;
    {
        int i=a->ref_end1,j=a->read_end1,state=0,nops=0;
        size_t maxops=(size_t)i+j+2;
        if(maxops>SIZE_MAX/sizeof(uint32_t)) {errno=ENOMEM;goto error;}
        uint32_t* ops=(uint32_t*)malloc(maxops*sizeof(uint32_t));
        if(!ops) {errno=ENOMEM;goto error;}
        while(i>=0&&j>=0) {
            uint8_t d=dir[(size_t)i*nq+j];
            int move=state?state:(d&3),op;
            if(!move) break;
            if(move==1) {op=0;--i;--j;state=0;}
            else if(move==2) {op=2;--i;state=(d>>2)&3;}
            else {op=1;--j;state=(d>>4)&3;}
            if(nops&&(ops[nops-1]&15)==(unsigned)op) ops[nops-1]+=16;
            else ops[nops++]=16|op;
        }
        a->ref_begin1=i+1;a->read_begin1=j+1;
        if(!(flag&7)||((flag&2)&&a->score1<filters)||
                ((flag&4)&&(a->ref_end1-a->ref_begin1>filterd||a->read_end1-a->read_begin1>filterd))) free(ops);
        else {
            for(int k=0;k<nops/2;++k) {uint32_t t=ops[k];ops[k]=ops[nops-1-k];ops[nops-1-k]=t;}
            a->cigar=ops;a->cigarLen=nops;
        }
    }
done:
    free(h);free(e);free(columns);free(dir);free(mstate);free(fstate);return a;
error:
    free(h);free(e);free(columns);free(dir);free(mstate);free(fstate);free(a->prefix_scores);free(a);return NULL;
}

static int ssw_avx2_valid_cigar(const s_profile* p,const int8_t* ref,const s_align* a,int go,int ge) {
    int q=a->read_begin1,r=a->ref_begin1;
    int64_t score=0;
    for(int i=0;i<a->cigarLen;++i) {
        unsigned n=a->cigar[i]>>4,op=a->cigar[i]&15;
        if(!n) return 0;
        if(op==0) {
            if(q<0||r<0||n>(unsigned)(a->read_end1+1-q)||n>(unsigned)(a->ref_end1+1-r)) return 0;
            for(unsigned j=0;j<n;++j) score+=p->mat[ref[r++]*p->n+p->read[q++]];
        } else {
            score-=go+(int64_t)(n-1)*ge;
            if(op==1) q+=n;else if(op==2) r+=n;else return 0;
        }
    }
    return score==a->score1&&q==a->read_end1+1&&r==a->ref_end1+1;
}
