from __future__ import division
import os, pysam
from collections import defaultdict
import numpy as np
import math

class Features:

    NAN = np.nan

    info_features_names = [ 'START_STOP_DIST', 'SVLEN', 'SVINSLEN', 'EDIT_DISTANCE',
                            'SV_REF_PREFIX_A_RATIO', 'SV_REF_PREFIX_C_RATIO', 'SV_REF_PREFIX_G_RATIO', 'SV_REF_PREFIX_T_RATIO', 'MAX_SV_REF_PREFIX_BASE_RATIO',
                            'SV_REF_SUFFIX_A_RATIO', 'SV_REF_SUFFIX_C_RATIO', 'SV_REF_SUFFIX_G_RATIO', 'SV_REF_SUFFIX_T_RATIO', 'MAX_SV_REF_SUFFIX_BASE_RATIO',
                            'LEFT_ANCHOR_A_RATIO', 'LEFT_ANCHOR_C_RATIO', 'LEFT_ANCHOR_G_RATIO', 'LEFT_ANCHOR_T_RATIO', 'MAX_LEFT_ANCHOR_BASE_RATIO',
                            'LEFT_FLANKING_A_RATIO_50', 'LEFT_FLANKING_C_RATIO_50', 'LEFT_FLANKING_G_RATIO_50', 'LEFT_FLANKING_T_RATIO_50', 'MAX_LEFT_FLANKING_BASE_RATIO_50',
                            'LEFT_FLANKING_A_RATIO_100', 'LEFT_FLANKING_C_RATIO_100', 'LEFT_FLANKING_G_RATIO_100', 'LEFT_FLANKING_T_RATIO_100', 'MAX_LEFT_FLANKING_BASE_RATIO_100',
                            'LEFT_FLANKING_A_RATIO_500', 'LEFT_FLANKING_C_RATIO_500', 'LEFT_FLANKING_G_RATIO_500', 'LEFT_FLANKING_T_RATIO_500', 'MAX_LEFT_FLANKING_BASE_RATIO_500',
                            'RIGHT_ANCHOR_A_RATIO', 'RIGHT_ANCHOR_C_RATIO', 'RIGHT_ANCHOR_G_RATIO', 'RIGHT_ANCHOR_T_RATIO', 'MAX_RIGHT_ANCHOR_BASE_RATIO',
                            'RIGHT_FLANKING_A_RATIO_50', 'RIGHT_FLANKING_C_RATIO_50', 'RIGHT_FLANKING_G_RATIO_50', 'RIGHT_FLANKING_T_RATIO_50', 'MAX_RIGHT_FLANKING_BASE_RATIO_50',
                            'RIGHT_FLANKING_A_RATIO_100', 'RIGHT_FLANKING_C_RATIO_100', 'RIGHT_FLANKING_G_RATIO_100', 'RIGHT_FLANKING_T_RATIO_100', 'MAX_RIGHT_FLANKING_BASE_RATIO_100',
                            'RIGHT_FLANKING_A_RATIO_500', 'RIGHT_FLANKING_C_RATIO_500', 'RIGHT_FLANKING_G_RATIO_500', 'RIGHT_FLANKING_T_RATIO_500', 'MAX_RIGHT_FLANKING_BASE_RATIO_500',
                            'INS_PREFIX_A_RATIO', 'INS_PREFIX_C_RATIO', 'INS_PREFIX_G_RATIO', 'INS_PREFIX_T_RATIO', 'MAX_INS_PREFIX_BASE_COUNT_RATIO',
                            'INS_SUFFIX_A_RATIO', 'INS_SUFFIX_C_RATIO', 'INS_SUFFIX_G_RATIO', 'INS_SUFFIX_T_RATIO', 'MAX_INS_SUFFIX_BASE_COUNT_RATIO',
                            'INS_SEQ_COV_PREFIX_LEN', 'INS_SEQ_COV_SUFFIX_LEN', 'EXP_ALT_READS_FREQ1', 'EXP_ALT_READS_FREQ2', 'HP_REF_LEN', 'HP_ALT_LEN' ]

    reads_features_names = ['AR1', 'AR1_ADJ', 'AR1C', 'AR1C_ADJ', 'AR1CmQ', 'AR1CMQ', 'AR1CHQ', 'AR1C_HQ_RATIO', 'AR1E', 'AR1E_RATIO',
                            'AR2', 'AR2_ADJ', 'AR2C', 'AR2C_ADJ', 'AR2CmQ', 'AR2CMQ', 'AR2CHQ', 'AR2C_HQ_RATIO', 'AR2E', 'AR2E_RATIO',
                            'AR1HPMODE', 'AR1CHPMODE', 'AR1CHPIQR', 'AR1HPMODE_AR1CHPMODE_DIFF', 'AR1HPMODE_ALTLEN_DIFF', 'AR1CHPMODE_ALTLEN_DIFF',
                            'AR1CHPmQ', 'AR1CHPMQ', 'AR1CHPAQ', 'AR1CHPSQ', 'AR1HP5PMR', 'AR1HP3PMR',
                            'MAXARCD', 'MAXARED',
                            'RR1', 'RR1C', 'RR1CmQ', 'RR1CMQ', 'RR1C_HQ_RATIO', 'RR1E', 'RR1E_RATIO',
                            'RR2', 'RR2C', 'RR2CmQ', 'RR2CMQ', 'RR2C_HQ_RATIO', 'RR2E', 'RR2E_RATIO', 'MAXRRCD', 'MAXRRED',
                            'OR1', 'OR2', 'OR1C', 'OR2C', 'OR1CHQ', 'OR2CHQ', 'OR1C_HQ_RATIO', 'OR2C_HQ_RATIO', 'OR1E', 'OR2E',
                            'NAR1', 'NAR2', 'NAR1C', 'NAR2C', 'NAR1CHQ', 'NAR2CHQ', 'NAR1C_HQ_RATIO', 'NAR2C_HQ_RATIO', 'NAR1E', 'NAR2E',
                            'ER',
                            'AR1CMSPAN_1', 'AR1CMSPAN_2', 'AR1CMHQSPAN_1', 'AR1CMHQSPAN_2',
                            'AR2CMSPAN_1', 'AR2CMSPAN_2', 'AR2CMHQSPAN_1', 'AR2CMHQSPAN_2',
                            'RR1CMSPAN_1', 'RR1CMSPAN_2', 'RR1CMHQSPAN_1', 'RR1CMHQSPAN_2',
                            'RR2CMSPAN_1', 'RR2CMSPAN_2', 'RR2CMHQSPAN_1', 'RR2CMHQSPAN_2',
                            'AR1_RR1_CAS_Z_SCORE', 'AR2_RR2_CAS_Z_SCORE', 
                            'AR1_OVER_RR1', 'AR2_OVER_RR2', 'AR1C_OVER_RR1C', 'AR2C_OVER_RR2C', 'AR1E_OVER_RR1E', 'AR2E_OVER_RR2E',
                            'AR1_OVER_OR1', 'AR2_OVER_OR2', 'AR1C_OVER_OR1C', 'AR2C_OVER_OR2C', 'AR1E_OVER_OR1E', 'AR2E_OVER_OR2E',
                            'AR1_OVER_NAR1', 'AR2_OVER_NAR2', 'AR1C_OVER_NAR1C', 'AR2C_OVER_NAR2C', 'AR1E_OVER_NAR1E', 'AR2E_OVER_NAR2E']

    fmt_features_names = [  'AXR1', 'AXR2', 'AXR1HQ', 'AXR2HQ',
                            'EXSS1_1', 'EXSS1_2', 'EXSS2_1', 'EXSS2_2',
                            'EXSS1_RATIO1', 'EXSS1_RATIO2', 'EXSS2_RATIO1', 'EXSS2_RATIO2',
                            'EXAS_EXRS_DIFF_TO_LEN',
                            'EXSSC1_IA_RATIO', 'EXSSC2_IA_RATIO', 'EXSSC1_IA_DIFF', 'EXSSC2_IA_DIFF',
                            'MEXL', 'mEXL', 'EXL',
                            'MDLF', 'MDSP', 'MDSF', 'MDRF', 'MDSP_OVER_MDLF', 'MDSF_OVER_MDRF',
                            'MDLFHQ', 'MDSPHQ', 'MDSFHQ', 'MDRFHQ', 'MDSP_OVER_MDLF_HQ', 'MDSF_OVER_MDRF_HQ',
                            'MDLC', 'MDRC', 'MDLCHQ', 'MDRCHQ',
                            'TD']

    stat_test_features_names = ['KS_PVAL', 'SIZE_NORM']

    dp_features_names = ['ASP1', 'ASP1HQ_1', 'ASP1HQ_2', 'ASP1HQ_1_RATIO', 'ASP1HQ_2_RATIO',
                         'ASP2', 'ASP2HQ_1', 'ASP2HQ_2', 'ASP2HQ_1_RATIO', 'ASP2HQ_2_RATIO',
                         'ASP1_ASP2_RATIO',
                         'ASP1mQ_1', 'ASP1mQ_2', 'ASP1MQ_1', 'ASP1MQ_2', 'ASP1SPAN_1', 'ASP1SPAN_2',
                         'ASP2mQ_1', 'ASP2mQ_2', 'ASP2MQ_1', 'ASP2MQ_2', 'ASP2SPAN_1', 'ASP2SPAN_2',
                         'ASP1_OVER_RSP1', 'ASP2_OVER_RSP2',
                         'ASP1_RSP1_1_NM_Z_SCORE', 'ASP1_RSP1_2_NM_Z_SCORE', 'ASP2_RSP2_1_NM_Z_SCORE', 'ASP2_RSP2_2_NM_Z_SCORE',
                         'RSP1', 'RSP1HQ_1', 'RSP1HQ_2',
                         'RSP2', 'RSP2HQ_1', 'RSP2HQ_2',
                         'RSP1mQ_1', 'RSP1mQ_2', 'RSP1MQ_1', 'RSP1MQ_2', 
                         'RSP2mQ_1', 'RSP2mQ_2', 'RSP2MQ_1', 'RSP2MQ_2',
                         'NSP1', 'NSP1HQ_1', 'NSP1HQ_2',
                         'NSP2', 'NSP2HQ_1', 'NSP2HQ_2',
                         'NSP1mQ_1', 'NSP1mQ_2', 'NSP1MQ_1', 'NSP1MQ_2',
                         'NSP2mQ_1', 'NSP2mQ_2', 'NSP2MQ_1', 'NSP2MQ_2',
                         'ASP1_NSP1_1_NM_Z_SCORE', 'ASP1_NSP1_2_NM_Z_SCORE', 'ASP2_NSP2_1_NM_Z_SCORE', 'ASP2_NSP2_2_NM_Z_SCORE',
                         'SSP1HQ_1', 'SSP1HQ_2',
                         'SSP2HQ_1', 'SSP2HQ_2',
                         'SSP1mQ_1', 'SSP1mQ_2', 'SSP1MQ_1', 'SSP1MQ_2',
                         'SSP2mQ_1', 'SSP2mQ_2', 'SSP2MQ_1', 'SSP2MQ_2',
                         'SSP1_RSP1_1_NM_Z_SCORE', 'SSP1_RSP1_2_NM_Z_SCORE', 'SSP2_RSP2_1_NM_Z_SCORE', 'SSP2_RSP2_2_NM_Z_SCORE'
    ]

    def get_feature_names(model_name):
        return Features.info_features_names + Features.fmt_features_names + Features.reads_features_names + \
            Features.stat_test_features_names + Features.dp_features_names

    def get_model_name(record, max_is, read_len):
        if Features.gt_as_homopolymer(record):
            return "HP"

        svtype_str = Features.get_svtype(record)
        sample = record.samples[0]

        if svtype_str == "DUP" and "INS_TO_DUP" in record.info:
            svtype_str = "INS_TO_DUP"
            if Features.get_svlen(record) > read_len-30:
                svtype_str += "_LARGE"

        if svtype_str == "DEL":
            if abs(Features.get_svlen(record)) >= max_is:
                svtype_str += "_LARGE"
            if 'EXL' not in sample:
                svtype_str += "_NOEXL"
        elif svtype_str == "DUP" and Features.get_svlen(record) > read_len-30:
            svtype_str += "_LARGE"
            if 'EXL' not in sample:
                svtype_str += "_NOEXL"

        return svtype_str

    def get_number_value(info, key, default, norm_factor = 1.0):
        if key in info:
            v = info[key]
        else:
            v = default
        if isinstance(v, list) or isinstance(v, tuple):
            return [float(x)/norm_factor for x in v]
        else:
            return float(v)/norm_factor

    def get_string_value(info, key, default):
        if key in info:
            return info[key]
        else:
            return default
        
    def get_string_list_value(info, key, default):
        if key in info:
            v = info[key]
            if isinstance(v, list) or isinstance(v, tuple):
                return [str(x) for x in v]
            else:
                return [str(v)]
        else:
            return default

    def generate_id(record):
        svinsseq = Features.get_svinsseq(record)
        aux_snps = Features.get_string_value(record.info, 'AUX_SNPS', "")
        aux_indels = Features.get_string_value(record.info, 'AUX_INDELS', "")
        return f"{record.chrom}:{record.pos}-{record.stop}:{Features.get_svtype(record)}:{Features.get_svlen(record)}:{hash(svinsseq)}:{hash(aux_snps)}:{hash(aux_indels)}"

    def get_svinsseq(record):
        if "<" not in record.alts[0]:
            return record.alts[0][1:]
        elif 'SVINSSEQ' in record.info:
            svinsseq = record.info['SVINSSEQ']
            if isinstance(svinsseq, list) or isinstance(svinsseq, tuple):
                svinsseq = svinsseq[0]
            return svinsseq
        elif "LEFT_SVINSSEQ" in record.info or "RIGHT_SVINSSEQ" in record.info:
            left_svinsseq = Features.get_string_value(record.info, 'LEFT_SVINSSEQ', "")
            right_svinsseq = Features.get_string_value(record.info, 'RIGHT_SVINSSEQ', "")
            if isinstance(left_svinsseq, list) or isinstance(left_svinsseq, tuple):
                left_svinsseq = left_svinsseq[0]
            if isinstance(right_svinsseq, list) or isinstance(right_svinsseq, tuple):
                right_svinsseq = right_svinsseq[0]
            return left_svinsseq + '-' + right_svinsseq
        return ""

    def get_svlen(record):
        svtype_str = Features.get_svtype(record)
        svinsseq = Features.get_svinsseq(record)
        if svtype_str in ["INS", "INS_TO_DUP", "DEL"]:
            svlen = len(svinsseq) - (record.stop - record.pos)
        elif svtype_str == "DUP":
            svlen = record.stop - record.pos + len(svinsseq)
        else:
            raise RuntimeError(f"Unexpected SVTYPE {svtype_str} for record {record.id}")

        if "AUX_INDELS" in record.info:
            aux_indels = Features.get_string_list_value(record.info, 'AUX_INDELS', [])
            for indel in aux_indels:
                sl = indel.split(':')
                svlen -= int(sl[1]) - int(sl[0]) # length of deletion
                svlen += len(sl[2]) # length of insertion
        return svlen

    def gt_as_homopolymer(record):
        return 'HP_GENOTYPED' in record.info

    def get_edit_distance(record, svinslen):
        edit_distance = record.stop - record.pos + svinslen
        if "AUX_SNPS" in record.info:
            aux_snps = Features.get_string_list_value(record.info, 'AUX_SNPS', [])
            edit_distance += len(aux_snps)
        if "AUX_INDELS" in record.info:
            aux_indels = Features.get_string_list_value(record.info, 'AUX_INDELS', [])
            for indel in aux_indels:
                sl = indel.split(':')
                edit_distance += int(sl[1]) - int(sl[0]) # length of deletion
                edit_distance += len(sl[2]) # length of insertion
        return edit_distance

    def get_svtype(record):
        if isinstance(record.info['SVTYPE'], list) or isinstance(record.info['SVTYPE'], tuple):
            return record.info['SVTYPE'][0]
        return record.info['SVTYPE']

    def skips_ml_genotyping(record):
        return Features.get_svtype(record).startswith('INV')

    def normalise(value, min, max):
        if max == min:
            return value - min
        return (value - min) / (max - min)
    
    def piecewise_normalise(value, minv, maxv):
        neg = value < 0
        value = abs(value)
        if value <= minv:
            ret_val = value/max(1, minv) * 0.25
        elif value > maxv:
            ret_val = value/max(1, maxv) * 0.25 + 0.75
        else:
            ret_val = 0.25 + (value - minv) / (maxv - minv) * 0.75
        if neg:
            return -ret_val
        return ret_val

    def calculate_z_score(mean1, stddev1, n1, mean2, stddev2, n2):
        if np.isnan(mean1) or np.isnan(mean2) or n1 == 0 or n2 == 0:
            return Features.NAN
        std_error = math.sqrt((stddev1**2 / n1) + (stddev2**2 / n2))
        if std_error == 0:
            std_error = 1
        z_score = (mean1 - mean2) / std_error
        return z_score

    def record_to_features(record, stats, feature_names = None):
        min_depth = get_stat(stats, 'min_depth', record.chrom)
        median_depth = get_stat(stats, 'median_depth', record.chrom)
        max_depth = get_stat(stats, 'max_depth', record.chrom)
        max_is = stats['max_is']['.']
        read_len = stats['read_len']['.']
        min_pairs_crossing_point = stats['min_pairs_crossing_gap']["0"]
        max_pairs_crossing_point = stats['max_pairs_crossing_gap']["0"]
        model_name = Features.get_model_name(record, max_is, read_len)

        features = dict()
        info = record.info
        sample = record.samples[0]
        svtype_str = Features.get_svtype(record)
        source_str = Features.get_string_value(info, 'SOURCE', "")
        features['START_STOP_DIST'] = record.stop - record.pos

        svlen = abs(Features.get_svlen(record))
        features['SVLEN'] = math.log1p(svlen)

        if feature_names is None:
            feature_names = Features.get_feature_names(model_name)

        svinsseq = Features.get_svinsseq(record)
        svinslen = len(svinsseq)
        features['SVINSLEN'] = svinslen

        edit_distance = Features.get_edit_distance(record, svinslen)
        features['EDIT_DISTANCE'] = edit_distance

        features['INS_SEQ_COV_PREFIX_LEN'] = 1
        features['INS_SEQ_COV_SUFFIX_LEN'] = 1
        if '-' in svinsseq:
            i = svinsseq.index('-')
            features['INS_SEQ_COV_PREFIX_LEN'] = i/len(svinsseq)
            features['INS_SEQ_COV_SUFFIX_LEN'] = (len(svinsseq)-i)/len(svinsseq)

        exp_alt_reads_freq1, exp_alt_reads_freq2 = info['EXP_ALT_READS_FREQ']
        features['EXP_ALT_READS_FREQ1'], features['EXP_ALT_READS_FREQ2'] = exp_alt_reads_freq1, exp_alt_reads_freq2

        if model_name == "HP":
            hp_ref_start, hp_ref_end = info['HP_REF_RANGE']
            features['HP_REF_LEN'] = hp_ref_end - hp_ref_start
            features['HP_ALT_LEN'] = features['HP_REF_LEN'] + Features.get_svlen(record)
        else:
             features['HP_REF_LEN'], features['HP_ALT_LEN'] = Features.NAN, Features.NAN

        left_anchor_base_count = info['LEFT_ANCHOR_BASE_COUNT'] 
        left_anchor_base_count_ratio = [x/max(1, sum(left_anchor_base_count)) for x in left_anchor_base_count]
        features['MAX_LEFT_ANCHOR_BASE_RATIO'] = max(left_anchor_base_count_ratio)
        features['LEFT_ANCHOR_A_RATIO'], features['LEFT_ANCHOR_C_RATIO'], features['LEFT_ANCHOR_G_RATIO'], features['LEFT_ANCHOR_T_RATIO'] = left_anchor_base_count_ratio

        left_flanking_base_count_50 = info['LEFT_FLANKING_BASE_COUNT_50']
        left_flanking_base_count_ratio_50 = [x/max(1, sum(left_flanking_base_count_50)) for x in left_flanking_base_count_50]
        features['MAX_LEFT_FLANKING_BASE_RATIO_50'] = max(left_flanking_base_count_ratio_50)
        features['LEFT_FLANKING_A_RATIO_50'], features['LEFT_FLANKING_C_RATIO_50'], features['LEFT_FLANKING_G_RATIO_50'], features['LEFT_FLANKING_T_RATIO_50'] = left_flanking_base_count_ratio_50

        left_flanking_base_count_100 = info['LEFT_FLANKING_BASE_COUNT_100']
        left_flanking_base_count_ratio_100 = [x/max(1, sum(left_flanking_base_count_100)) for x in left_flanking_base_count_100]
        features['MAX_LEFT_FLANKING_BASE_RATIO_100'] = max(left_flanking_base_count_ratio_100)
        features['LEFT_FLANKING_A_RATIO_100'], features['LEFT_FLANKING_C_RATIO_100'], features['LEFT_FLANKING_G_RATIO_100'], features['LEFT_FLANKING_T_RATIO_100'] = left_flanking_base_count_ratio_100

        left_flanking_base_count_500 = info['LEFT_FLANKING_BASE_COUNT_500']
        left_flanking_base_count_ratio_500 = [x/max(1, sum(left_flanking_base_count_500)) for x in left_flanking_base_count_500]
        features['MAX_LEFT_FLANKING_BASE_RATIO_500'] = max(left_flanking_base_count_ratio_500)
        features['LEFT_FLANKING_A_RATIO_500'], features['LEFT_FLANKING_C_RATIO_500'], features['LEFT_FLANKING_G_RATIO_500'], features['LEFT_FLANKING_T_RATIO_500'] = left_flanking_base_count_ratio_500

        right_anchor_base_count = info['RIGHT_ANCHOR_BASE_COUNT']
        right_anchor_base_count_ratio = [x/max(1, sum(right_anchor_base_count)) for x in right_anchor_base_count]
        features['MAX_RIGHT_ANCHOR_BASE_RATIO'] = max(right_anchor_base_count_ratio)
        features['RIGHT_ANCHOR_A_RATIO'], features['RIGHT_ANCHOR_C_RATIO'], features['RIGHT_ANCHOR_G_RATIO'], features['RIGHT_ANCHOR_T_RATIO'] = right_anchor_base_count_ratio

        right_flanking_base_count_50 = info['RIGHT_FLANKING_BASE_COUNT_50']
        right_flanking_base_count_ratio_50 = [x/max(1, sum(right_flanking_base_count_50)) for x in right_flanking_base_count_50]
        features['MAX_RIGHT_FLANKING_BASE_RATIO_50'] = max(right_flanking_base_count_ratio_50)
        features['RIGHT_FLANKING_A_RATIO_50'], features['RIGHT_FLANKING_C_RATIO_50'], features['RIGHT_FLANKING_G_RATIO_50'], features['RIGHT_FLANKING_T_RATIO_50'] = right_flanking_base_count_ratio_50

        right_flanking_base_count_100 = info['RIGHT_FLANKING_BASE_COUNT_100']
        right_flanking_base_count_ratio_100 = [x/max(1, sum(right_flanking_base_count_100)) for x in right_flanking_base_count_100]
        features['MAX_RIGHT_FLANKING_BASE_RATIO_100'] = max(right_flanking_base_count_ratio_100)
        features['RIGHT_FLANKING_A_RATIO_100'], features['RIGHT_FLANKING_C_RATIO_100'], features['RIGHT_FLANKING_G_RATIO_100'], features['RIGHT_FLANKING_T_RATIO_100'] = right_flanking_base_count_ratio_100

        right_flanking_base_count_500 = info['RIGHT_FLANKING_BASE_COUNT_500']
        right_flanking_base_count_ratio_500 = [x/max(1, sum(right_flanking_base_count_500)) for x in right_flanking_base_count_500]
        features['MAX_RIGHT_FLANKING_BASE_RATIO_500'] = max(right_flanking_base_count_ratio_500)
        features['RIGHT_FLANKING_A_RATIO_500'], features['RIGHT_FLANKING_C_RATIO_500'], features['RIGHT_FLANKING_G_RATIO_500'], features['RIGHT_FLANKING_T_RATIO_500'] = right_flanking_base_count_ratio_500

        sv_ref_prefix_base_count = info['SV_REF_PREFIX_BASE_COUNT']
        sv_ref_prefix_base_count_ratio = [x/max(1, sum(sv_ref_prefix_base_count)) for x in sv_ref_prefix_base_count]
        features['MAX_SV_REF_PREFIX_BASE_RATIO'] = max(sv_ref_prefix_base_count_ratio)
        features['SV_REF_PREFIX_A_RATIO'], features['SV_REF_PREFIX_C_RATIO'], features['SV_REF_PREFIX_G_RATIO'], features['SV_REF_PREFIX_T_RATIO'] = sv_ref_prefix_base_count_ratio

        sv_ref_suffix_base_count = info['SV_REF_SUFFIX_BASE_COUNT']
        sv_ref_suffix_base_count_ratio = [x/max(1, sum(sv_ref_suffix_base_count)) for x in sv_ref_suffix_base_count]
        features['MAX_SV_REF_SUFFIX_BASE_RATIO'] = max(sv_ref_suffix_base_count_ratio)
        features['SV_REF_SUFFIX_A_RATIO'], features['SV_REF_SUFFIX_C_RATIO'], features['SV_REF_SUFFIX_G_RATIO'], features['SV_REF_SUFFIX_T_RATIO'] = sv_ref_suffix_base_count_ratio

        ins_prefix_base_count = info['INS_PREFIX_BASE_COUNT']
        ins_prefix_base_count_ratio = [x/max(1, sum(ins_prefix_base_count)) for x in ins_prefix_base_count]
        features['MAX_INS_PREFIX_BASE_COUNT_RATIO'] = max(ins_prefix_base_count_ratio)
        features['INS_PREFIX_A_RATIO'], features['INS_PREFIX_C_RATIO'], features['INS_PREFIX_G_RATIO'], features['INS_PREFIX_T_RATIO'] = ins_prefix_base_count_ratio

        ins_suffix_base_count = info['INS_SUFFIX_BASE_COUNT']
        ins_suffix_base_count_ratio = [x/max(1, sum(ins_suffix_base_count)) for x in ins_suffix_base_count]
        features['MAX_INS_SUFFIX_BASE_COUNT_RATIO'] = max(ins_suffix_base_count_ratio)
        features['INS_SUFFIX_A_RATIO'], features['INS_SUFFIX_C_RATIO'], features['INS_SUFFIX_G_RATIO'], features['INS_SUFFIX_T_RATIO'] = ins_suffix_base_count_ratio

        features['TD'] = Features.get_number_value(sample, 'TD', 0)

        ar1 = Features.get_number_value(sample, 'AR1', 0)
        ar1c = Features.get_number_value(sample, 'AR1C', 0)
        arc1hq = Features.get_number_value(sample, 'AR1CHQ', 0)
        ar1e = Features.get_number_value(sample, 'AR1E', 0)
        ar1_adj = ar1
        ar1c_adj = ar1c
        if exp_alt_reads_freq1 > 0:
            ar1_adj = ar1/exp_alt_reads_freq1
            ar1c_adj = ar1c/exp_alt_reads_freq1
        ar1cas = Features.get_number_value(sample, 'AR1CAS', Features.NAN)
        ar1css = Features.get_number_value(sample, 'AR1CSS', Features.NAN)
        features['AR1'] = Features.piecewise_normalise(ar1, min_depth, max_depth)
        features['AR1C'] = Features.piecewise_normalise(ar1c, min_depth, max_depth)
        features['AR1_ADJ'] = Features.piecewise_normalise(ar1_adj, min_depth, max_depth)
        features['AR1C_ADJ'] = Features.piecewise_normalise(ar1c_adj, min_depth, max_depth)
        features['AR1CmQ'] = Features.get_number_value(sample, 'AR1CmQ', Features.NAN)
        features['AR1CMQ'] = Features.get_number_value(sample, 'AR1CMQ', Features.NAN)
        features['AR1CHQ'] = Features.piecewise_normalise(arc1hq, min_depth, max_depth)
        features['AR1C_HQ_RATIO'] = arc1hq/max(1, ar1c)
        features['AR1E'] = Features.piecewise_normalise(ar1e, min_depth, max_depth)
        features['AR1E_RATIO'] = ar1e/max(1, ar1c)
        features['AR1CMSPAN_1'], features['AR1CMSPAN_2'] = Features.get_number_value(sample, 'AR1CMSPAN', [0, 0], max_is)
        features['AR1CMHQSPAN_1'], features['AR1CMHQSPAN_2'] = Features.get_number_value(sample, 'AR1CMHQSPAN', [0, 0], max_is)

        features['AR1HPMODE'] = Features.get_number_value(sample, 'AR1HPMODE', Features.NAN)
        features['AR1CHPMODE'] = Features.get_number_value(sample, 'AR1CHPMODE', Features.NAN)
        features['AR1HPMODE_AR1CHPMODE_DIFF'] = features['AR1HPMODE'] - features['AR1CHPMODE']
        features['AR1HPMODE_ALTLEN_DIFF'] = features['AR1HPMODE'] - features['HP_ALT_LEN']
        features['AR1CHPMODE_ALTLEN_DIFF'] = features['AR1CHPMODE'] - features['HP_ALT_LEN']
        features['AR1CHPIQR'] = Features.get_number_value(sample, 'AR1CHPIQR', Features.NAN)
        features['AR1CHPmQ'] = Features.get_number_value(sample, 'AR1CHPmQ', Features.NAN)
        features['AR1CHPMQ'] = Features.get_number_value(sample, 'AR1CHPMQ', Features.NAN)
        features['AR1CHPAQ'] = Features.get_number_value(sample, 'AR1CHPAQ', Features.NAN)
        features['AR1CHPSQ'] = Features.get_number_value(sample, 'AR1CHPSQ', Features.NAN)
        features['AR1HP5PMR'] = Features.get_number_value(sample, 'AR1HP5PMR', Features.NAN)
        features['AR1HP3PMR'] = Features.get_number_value(sample, 'AR1HP3PMR', Features.NAN)

        ar2 = Features.get_number_value(sample, 'AR2', 0)
        ar2c = Features.get_number_value(sample, 'AR2C', 0)
        ar2_adj = ar2
        ar2c_adj = ar2c
        if exp_alt_reads_freq2 > 0:
            ar2_adj = ar2/exp_alt_reads_freq2
            ar2c_adj = ar2c/exp_alt_reads_freq2
        ar2cas = Features.get_number_value(sample, 'AR2CAS', Features.NAN)
        ar2css = Features.get_number_value(sample, 'AR2CSS', Features.NAN)
        features['AR2'] = Features.piecewise_normalise(ar2, min_depth, max_depth)
        features['AR2C'] = Features.piecewise_normalise(ar2c, min_depth, max_depth)
        features['AR2_ADJ'] = Features.piecewise_normalise(ar2_adj, min_depth, max_depth)
        features['AR2C_ADJ'] = Features.piecewise_normalise(ar2c_adj, min_depth, max_depth)
        features['AR2CmQ'] = Features.get_number_value(sample, 'AR2CmQ', Features.NAN)
        features['AR2CMQ'] = Features.get_number_value(sample, 'AR2CMQ', Features.NAN)
        arc2hq = Features.get_number_value(sample, 'AR2CHQ', 0)
        features['AR2CHQ'] = Features.piecewise_normalise(arc2hq, min_depth, max_depth)
        features['AR2C_HQ_RATIO'] = arc2hq/max(1, ar2c)
        ar2e = Features.get_number_value(sample, 'AR2E', 0)
        features['AR2E'] = Features.piecewise_normalise(ar2e, min_depth, max_depth)
        features['AR2E_RATIO'] = ar2e/max(1, ar2c)
        features['AR2CMSPAN_1'], features['AR2CMSPAN_2'] = Features.get_number_value(sample, 'AR2CMSPAN', [0, 0], max_is)
        features['AR2CMHQSPAN_1'], features['AR2CMHQSPAN_2'] = Features.get_number_value(sample, 'AR2CMHQSPAN', [0, 0], max_is)

        ar1cf = Features.get_number_value(sample, 'AR1CF', 0, max(1, ar1c))
        ar1cr = Features.get_number_value(sample, 'AR1CR', 0, max(1, ar1c))
        ar2cf = Features.get_number_value(sample, 'AR2CF', 0, max(1, ar2c))
        ar2cr = Features.get_number_value(sample, 'AR2CR', 0, max(1, ar2c))
        features['ARCF'] = ar1cf + ar2cf
        features['ARCR'] = ar1cr + ar2cr
        features['MAXARCD'] = max(features['ARCF'], features['ARCR'])

        ar1ef = Features.get_number_value(sample, 'AR1EF', 0, max(1, ar1e))
        ar1er = Features.get_number_value(sample, 'AR1ER', 0, max(1, ar1e))
        ar2ef = Features.get_number_value(sample, 'AR2EF', 0, max(1, ar2e))
        ar2er = Features.get_number_value(sample, 'AR2ER', 0, max(1, ar2e))
        features['AREF'] = ar1ef + ar2ef
        features['ARER'] = ar1er + ar2er
        features['MAXARED'] = max(features['AREF'], features['ARER'])

        or1 = Features.get_number_value(sample, 'OR1', 0)
        or1c = Features.get_number_value(sample, 'OR1C', 0)
        or1chq = Features.get_number_value(sample, 'OR1CHQ', 0)
        or1e = Features.get_number_value(sample, 'OR1E', 0)
        or2 = Features.get_number_value(sample, 'OR2', 0)
        or2c = Features.get_number_value(sample, 'OR2C', 0)
        or2chq = Features.get_number_value(sample, 'OR2CHQ', 0)
        or2e = Features.get_number_value(sample, 'OR2E', 0)
        features['OR1'] = Features.piecewise_normalise(or1, min_depth, max_depth)
        features['OR1C'] = Features.piecewise_normalise(or1c, min_depth, max_depth)
        features['OR1CHQ'] = Features.piecewise_normalise(or1chq, min_depth, max_depth)
        features['OR1C_HQ_RATIO'] = or1chq/max(1, or1c)
        features['OR1E'] = Features.piecewise_normalise(or1e, min_depth, max_depth)
        features['OR2'] = Features.piecewise_normalise(or2, min_depth, max_depth)
        features['OR2C'] = Features.piecewise_normalise(or2c, min_depth, max_depth)
        features['OR2CHQ'] = Features.piecewise_normalise(or2chq, min_depth, max_depth)
        features['OR2C_HQ_RATIO'] = or2chq/max(1, or2c)
        features['OR2E'] = Features.piecewise_normalise(or2e, min_depth, max_depth)

        rr1 = Features.get_number_value(sample, 'RR1', 0)
        rr1c = Features.get_number_value(sample, 'RR1C', 0)
        rr1chq = Features.get_number_value(sample, 'RR1CHQ', 0)
        rr1e = Features.get_number_value(sample, 'RR1E', 0)
        rr2 = Features.get_number_value(sample, 'RR2', 0)
        rr2c = Features.get_number_value(sample, 'RR2C', 0)
        rr2chq = Features.get_number_value(sample, 'RR2CHQ', 0)
        rr2e = Features.get_number_value(sample, 'RR2E', 0)

        rr1cas = Features.get_number_value(sample, 'RR1CAS', Features.NAN)
        rr1css = Features.get_number_value(sample, 'RR1CSS', Features.NAN)
        features['RR1'] = Features.piecewise_normalise(rr1, min_depth, max_depth)
        features['RR1C'] = Features.piecewise_normalise(rr1c, min_depth, max_depth)
        features['RR1CmQ'] = Features.get_number_value(sample, 'RR1CmQ', Features.NAN)
        features['RR1CMQ'] = Features.get_number_value(sample, 'RR1CMQ', Features.NAN)
        features['RR1C_HQ_RATIO'] = rr1chq/max(1, rr1c)
        features['RR1E'] = Features.piecewise_normalise(rr1e, min_depth, max_depth)
        features['RR1E_RATIO'] = rr1e/max(1, rr1c)
        features['RR1CMSPAN_1'], features['RR1CMSPAN_2'] = Features.get_number_value(sample, 'RR1CMSPAN', [0, 0], max_is)
        features['RR1CMHQSPAN_1'], features['RR1CMHQSPAN_2'] = Features.get_number_value(sample, 'RR1CMHQSPAN', [0, 0], max_is)

        rr2cas = Features.get_number_value(sample, 'RR2CAS', Features.NAN)
        rr2css = Features.get_number_value(sample, 'RR2CSS', Features.NAN)
        features['RR2'] = Features.piecewise_normalise(rr2, min_depth, max_depth)
        features['RR2C'] = Features.piecewise_normalise(rr2c, min_depth, max_depth)
        features['RR2CmQ'] = Features.get_number_value(sample, 'RR2CmQ', Features.NAN)
        features['RR2CMQ'] = Features.get_number_value(sample, 'RR2CMQ', Features.NAN)
        features['RR2C_HQ_RATIO'] = rr2chq/max(1, rr2c)
        features['RR2E'] = Features.piecewise_normalise(rr2e, min_depth, max_depth)
        features['RR2E_RATIO'] = rr2e/max(1, rr2c)
        features['RR2CMSPAN_1'], features['RR2CMSPAN_2'] = Features.get_number_value(sample, 'RR2CMSPAN', [0, 0], max_is)
        features['RR2CMHQSPAN_1'], features['RR2CMHQSPAN_2'] = Features.get_number_value(sample, 'RR2CMHQSPAN', [0, 0], max_is)

        rr1cf = Features.get_number_value(sample, 'RR1CF', 0)
        rr1cr = Features.get_number_value(sample, 'RR1CR', 0)
        rr2cf = Features.get_number_value(sample, 'RR2CF', 0)
        rr2cr = Features.get_number_value(sample, 'RR2CR', 0)
        rr1cf_ratio = rr1cf/max(1, rr1cf + rr1cr)
        rr1cr_ratio = rr1cr/max(1, rr1cf + rr1cr)
        rr2cf_ratio = rr2cf/max(1, rr2cf + rr2cr)
        rr2cr_ratio = rr2cr/max(1, rr2cf + rr2cr)
        features['RRCF'] = rr1cf_ratio + rr2cf_ratio
        features['RRCR'] = rr1cr_ratio + rr2cr_ratio
        features['MAXRRCD'] = max(features['RRCF'], features['RRCR'])

        rr1ef = Features.get_number_value(sample, 'RR1EF', 0, max(1, rr1e))
        rr1er = Features.get_number_value(sample, 'RR1ER', 0, max(1, rr1e))
        rr2ef = Features.get_number_value(sample, 'RR2EF', 0, max(1, rr2e))
        rr2er = Features.get_number_value(sample, 'RR2ER', 0, max(1, rr2e))
        features['RREF'] = rr1ef + rr2ef
        features['RRER'] = rr1er + rr2er
        features['MAXRRED'] = max(features['RREF'], features['RRER'])

        nar1 = rr1 + or1
        nar1c = rr1c + or1c
        nar1chq = rr1chq + or1chq
        nar1e = rr1e + or1e
        nar2 = rr2 + or2
        nar2c = rr2c + or2c
        nar2chq = rr2chq + or2chq
        nar2e = rr2e + or2e
        features['NAR1'] = Features.piecewise_normalise(nar1, min_depth, max_depth)
        features['NAR1C'] = Features.piecewise_normalise(nar1c, min_depth, max_depth)
        features['NAR1CHQ'] = Features.piecewise_normalise(nar1chq, min_depth, max_depth)
        features['NAR1C_HQ_RATIO'] = nar1chq/max(1, nar1c)
        features['NAR1E'] = Features.piecewise_normalise(nar1e, min_depth, max_depth)
        features['NAR2'] = Features.piecewise_normalise(nar2, min_depth, max_depth)
        features['NAR2C'] = Features.piecewise_normalise(nar2c, min_depth, max_depth)
        features['NAR2CHQ'] = Features.piecewise_normalise(nar2chq, min_depth, max_depth)
        features['NAR2C_HQ_RATIO'] = nar2chq/max(1, nar2c)
        features['NAR2E'] = Features.piecewise_normalise(nar2e, min_depth, max_depth)

        er = Features.get_number_value(sample, 'ER', 0)
        features['ER'] = Features.piecewise_normalise(er, min_depth, max_depth)

        features['AR1_RR1_CAS_Z_SCORE'] = Features.calculate_z_score(ar1cas, ar1css, ar1c, rr1cas, rr1css, rr1c)

        if 'AR2' not in sample:
            ar2 = ar1
            ar2c = ar1c
            ar2cas = ar1cas
            ar2css = ar1css
            ar2e = ar1e
        if 'RR2' not in sample:
            rr2 = rr1
            rr2c = rr1c
            rr2cas = rr1cas
            rr2css = rr1css
            rr2e = rr1e

        features['AR2_RR2_CAS_Z_SCORE'] = Features.calculate_z_score(ar2cas, ar2css, ar2c, rr2cas, rr2css, rr2c)

        features['AR1_OVER_RR1'] = ar1/max(1, ar1+rr1)
        features['AR2_OVER_RR2'] = ar2/max(1, ar2+rr2)
        features['AR1C_OVER_RR1C'] = ar1c/max(1, ar1c+rr1c)
        features['AR2C_OVER_RR2C'] = ar2c/max(1, ar2c+rr2c)
        features['AR1E_OVER_RR1E'] = ar1e/max(1, ar1e+rr1e)
        features['AR2E_OVER_RR2E'] = ar2e/max(1, ar2e+rr2e)
        features['AR1_OVER_OR1'] = ar1/max(1, ar1+or1)
        features['AR2_OVER_OR2'] = ar2/max(1, ar2+or2)
        features['AR1C_OVER_OR1C'] = ar1c/max(1, ar1c+or1c)
        features['AR2C_OVER_OR2C'] = ar2c/max(1, ar2c+or2c)
        features['AR1E_OVER_OR1E'] = ar1e/max(1, ar1e+or1e)
        features['AR2E_OVER_OR2E'] = ar2e/max(1, ar2e+or2e)
        features['AR1_OVER_NAR1'] = ar1/max(1, ar1+nar1)
        features['AR2_OVER_NAR2'] = ar2/max(1, ar2+nar2)
        features['AR1C_OVER_NAR1C'] = ar1c/max(1, ar1c+nar1c)
        features['AR2C_OVER_NAR2C'] = ar2c/max(1, ar2c+nar2c)
        features['AR1E_OVER_NAR1E'] = ar1e/max(1, ar1e+nar1e)
        features['AR2E_OVER_NAR2E'] = ar2e/max(1, ar2e+nar2e)

        md = Features.get_number_value(sample, 'MD', [0, 0, 0, 0])
        features['MDLF'] = Features.piecewise_normalise(md[0], min_depth, max_depth)
        features['MDSP'] = Features.piecewise_normalise(md[1], min_depth, max_depth)
        features['MDSF'] = Features.piecewise_normalise(md[2], min_depth, max_depth)
        features['MDRF'] = Features.piecewise_normalise(md[3], min_depth, max_depth)
        features['MDSP_OVER_MDLF'] = md[1]/max(1, md[0])
        features['MDSF_OVER_MDRF'] = md[2]/max(1, md[3])

        mdhq = Features.get_number_value(sample, 'MDHQ', [0, 0, 0, 0])
        features['MDLFHQ'] = Features.piecewise_normalise(mdhq[0], min_depth, max_depth)
        features['MDSPHQ'] = Features.piecewise_normalise(mdhq[1], min_depth, max_depth)
        features['MDSFHQ'] = Features.piecewise_normalise(mdhq[2], min_depth, max_depth)
        features['MDRFHQ'] = Features.piecewise_normalise(mdhq[3], min_depth, max_depth)
        features['MDSP_OVER_MDLF_HQ'] = Features.piecewise_normalise(mdhq[1]-mdhq[0], min_depth, max_depth)
        features['MDSF_OVER_MDRF_HQ'] = Features.piecewise_normalise(mdhq[2]-mdhq[3], min_depth, max_depth)

        clmd = Features.get_number_value(sample, 'CLMD', [0, 0])
        features['MDLC'] = Features.piecewise_normalise(clmd[0], min_depth, max_depth)
        features['MDRC'] = Features.piecewise_normalise(clmd[1], min_depth, max_depth)

        clmdhq = Features.get_number_value(sample, 'CLMDHQ', [0, 0])
        features['MDLCHQ'] = Features.piecewise_normalise(clmdhq[0], min_depth, max_depth)
        features['MDRCHQ'] = Features.piecewise_normalise(clmdhq[1], min_depth, max_depth)

        features['KS_PVAL'] = Features.get_number_value(sample, 'KSPVAL', Features.NAN)
        features['SIZE_NORM'] = Features.NAN
        if 'MAXSIZE' in sample:
            min_size = float(sample['MINSIZE'])
            max_size = float(sample['MAXSIZE'])
            features['SIZE_NORM'] = Features.normalise(svlen/2, min_size, max_size)

        if svtype_str == "DEL":
            min_is_to_become_disc = int(max(0, max_is-svlen))
            min_disc_pairs = stats['min_pairs_crossing_gap'][str(min_is_to_become_disc)]
            max_disc_pairs = stats['max_pairs_crossing_gap'][str(min_is_to_become_disc)]
        elif svtype_str == "INS" and source_str in ("DE_NOVO_ASSEMBLY", "REFERENCE_GUIDED_ASSEMBLY"):
            min_inslen = int(min(max_is, svinslen))
            min_disc_pairs = stats['min_disc_pairs_by_insertion_size'][str(min_inslen)]
            max_disc_pairs = stats['max_disc_pairs_by_insertion_size'][str(min_inslen)]
        else:
            min_disc_pairs = min_pairs_crossing_point
            max_disc_pairs = max_pairs_crossing_point

        asp1 = Features.get_number_value(sample, 'ASP1', 0)
        asp1hq_1, asp1hq_2 = Features.get_number_value(sample, 'ASP1HQ', [0, 0])
        asp1nma_1, asp1nma_2 = Features.get_number_value(sample, 'ASP1NMA', [Features.NAN, Features.NAN])
        asp1nms_1, asp1nms_2 = Features.get_number_value(sample, 'ASP1NMS', [Features.NAN, Features.NAN])
        features['ASP1'] = Features.piecewise_normalise(asp1, min_disc_pairs, max_disc_pairs)
        features['ASP1HQ_1'] = Features.piecewise_normalise(asp1hq_1, min_disc_pairs, max_disc_pairs)
        features['ASP1HQ_2'] = Features.piecewise_normalise(asp1hq_2, min_disc_pairs, max_disc_pairs)
        features['ASP1HQ_1_RATIO'], features['ASP1HQ_2_RATIO'] = asp1hq_1/max(1, asp1), asp1hq_2/max(1, asp1)
        features['ASP1mQ_1'], features['ASP1mQ_2'] = Features.get_number_value(sample, 'ASP1mQ', [Features.NAN, Features.NAN])
        features['ASP1MQ_1'], features['ASP1MQ_2'] = Features.get_number_value(sample, 'ASP1MQ', [Features.NAN, Features.NAN])

        asp2 = Features.get_number_value(sample, 'ASP2', 0)
        asp2hq_1, asp2hq_2 = Features.get_number_value(sample, 'ASP2HQ', [0, 0])
        asp2nma_1, asp2nma_2 = Features.get_number_value(sample, 'ASP2NMA', [Features.NAN, Features.NAN])
        asp2nms_1, asp2nms_2 = Features.get_number_value(sample, 'ASP2NMS', [Features.NAN, Features.NAN])
        features['ASP2'] = Features.piecewise_normalise(asp2, min_disc_pairs, max_disc_pairs)
        features['ASP2HQ_1'] = Features.piecewise_normalise(asp2hq_1, min_disc_pairs, max_disc_pairs)
        features['ASP2HQ_2'] = Features.piecewise_normalise(asp2hq_2, min_disc_pairs, max_disc_pairs)
        features['ASP2HQ_1_RATIO'], features['ASP2HQ_2_RATIO'] = asp2hq_1/max(1, asp2), asp2hq_2/max(1, asp2)
        features['ASP2mQ_1'], features['ASP2mQ_2'] = Features.get_number_value(sample, 'ASP2mQ', [Features.NAN, Features.NAN])
        features['ASP2MQ_1'], features['ASP2MQ_2'] = Features.get_number_value(sample, 'ASP2MQ', [Features.NAN, Features.NAN])

        features['ASP1_ASP2_RATIO'] = max(asp1, asp2)/max(1, asp1+asp2)

        asp1span_1, asp1span_2 = Features.get_number_value(sample, 'ASP1SPAN', [0, 0])
        asp2span_1, asp2span_2 = Features.get_number_value(sample, 'ASP2SPAN', [0, 0])
        features['ASP1SPAN_1'], features['ASP2SPAN_2'] = asp1span_1/max_is, asp2span_2/max_is
        if svtype_str == "INS":
            features['ASP1SPAN_2'] = asp1span_2/max(1, max_is, svinslen)
            features['ASP2SPAN_1'] = asp2span_1/max(1, max_is, svinslen)
        else:
            features['ASP1SPAN_2'] = asp1span_2/max_is
            features['ASP2SPAN_1'] = asp2span_1/max_is

        rsp1 = Features.get_number_value(sample, 'RSP1', 0)
        rsp1hq_1, rsp1hq_2 = Features.get_number_value(sample, 'RSP1HQ', [0, 0])
        rsp1nma_1, rsp1nma_2 = Features.get_number_value(sample, 'RSP1NMA', [Features.NAN, Features.NAN])
        rsp1nms_1, rsp1nms_2 = Features.get_number_value(sample, 'RSP1NMS', [Features.NAN, Features.NAN])
        features['RSP1'] = Features.piecewise_normalise(rsp1, min_disc_pairs, max_disc_pairs)
        features['RSP1HQ_1']= Features.piecewise_normalise(rsp1hq_1, min_disc_pairs, max_disc_pairs)
        features['RSP1HQ_2'] = Features.piecewise_normalise(rsp1hq_2, min_disc_pairs, max_disc_pairs)
        features['RSP1HQ_1_RATIO'], features['RSP1HQ_2_RATIO'] = rsp1hq_1/max(1, rsp1), rsp1hq_2/max(1, rsp1)
        features['RSP1mQ_1'], features['RSP1mQ_2'] = Features.get_number_value(sample, 'RSP1mQ', [Features.NAN, Features.NAN])
        features['RSP1MQ_1'], features['RSP1MQ_2'] = Features.get_number_value(sample, 'RSP1MQ', [Features.NAN, Features.NAN])

        rsp2 = Features.get_number_value(sample, 'RSP2', 0)
        rsp2hq_1, rsp2hq_2 = Features.get_number_value(sample, 'RSP2HQ', [0, 0])
        rsp2nma_1, rsp2nma_2 = Features.get_number_value(sample, 'RSP2NMA', [Features.NAN, Features.NAN])
        rsp2nms_1, rsp2nms_2 = Features.get_number_value(sample, 'RSP2NMS', [Features.NAN, Features.NAN])
        features['RSP2'] = Features.piecewise_normalise(rsp2, min_disc_pairs, max_disc_pairs)
        features['RSP2HQ_1'] = Features.piecewise_normalise(rsp2hq_1, min_disc_pairs, max_disc_pairs)
        features['RSP2HQ_2'] = Features.piecewise_normalise(rsp2hq_2, min_disc_pairs, max_disc_pairs)
        features['RSP2HQ_1_RATIO'], features['RSP2HQ_2_RATIO'] = rsp2hq_1/max(1, rsp2), rsp2hq_2/max(1, rsp2)
        features['RSP2mQ_1'], features['RSP2mQ_2'] = Features.get_number_value(sample, 'RSP2mQ', [Features.NAN, Features.NAN])
        features['RSP2MQ_1'], features['RSP2MQ_2'] = Features.get_number_value(sample, 'RSP2MQ', [Features.NAN, Features.NAN])

        nsp1 = Features.get_number_value(sample, 'NSP1', 0)
        nsp1hq_1, nsp1hq_2 = Features.get_number_value(sample, 'NSP1HQ', [0, 0])
        nsp1nma_1, nsp1nma_2 = Features.get_number_value(sample, 'NSP1NMA', [Features.NAN, Features.NAN])
        nsp1nms_1, nsp1nms_2 = Features.get_number_value(sample, 'NSP1NMS', [Features.NAN, Features.NAN])
        features['NSP1'] = Features.piecewise_normalise(nsp1, min_disc_pairs, max_disc_pairs)
        features['NSP1HQ_1'] = Features.piecewise_normalise(nsp1hq_1, min_disc_pairs, max_disc_pairs)
        features['NSP1HQ_2'] = Features.piecewise_normalise(nsp1hq_2, min_disc_pairs, max_disc_pairs)
        features['NSP1HQ_1_RATIO'], features['NSP1HQ_2_RATIO'] = nsp1hq_1/max(1, nsp1), nsp1hq_2/max(1, nsp1)
        features['NSP1mQ_1'], features['NSP1mQ_2'] = Features.get_number_value(sample, 'NSP1mQ', [Features.NAN, Features.NAN])
        features['NSP1MQ_1'], features['NSP1MQ_2'] = Features.get_number_value(sample, 'NSP1MQ', [Features.NAN, Features.NAN])

        nsp2 = Features.get_number_value(sample, 'NSP2', 0)
        nsp2hq_1, nsp2hq_2 = Features.get_number_value(sample, 'NSP2HQ', [0, 0])
        nsp2nma_1, nsp2nma_2 = Features.get_number_value(sample, 'NSP2NMA', [Features.NAN, Features.NAN])
        nsp2nms_1, nsp2nms_2 = Features.get_number_value(sample, 'NSP2NMS', [Features.NAN, Features.NAN])
        features['NSP2'] = Features.piecewise_normalise(nsp2, min_disc_pairs, max_disc_pairs)
        features['NSP2HQ_1'] = Features.piecewise_normalise(nsp2hq_1, min_disc_pairs, max_disc_pairs)
        features['NSP2HQ_2'] = Features.piecewise_normalise(nsp2hq_2, min_disc_pairs, max_disc_pairs)
        features['NSP2HQ_1_RATIO'], features['NSP2HQ_2_RATIO'] = nsp2hq_1/max(1, nsp2), nsp2hq_2/max(1, nsp2)
        features['NSP2mQ_1'], features['NSP2mQ_2'] = Features.get_number_value(sample, 'NSP2mQ', [Features.NAN, Features.NAN])
        features['NSP2MQ_1'], features['NSP2MQ_2'] = Features.get_number_value(sample, 'NSP2MQ', [Features.NAN, Features.NAN])

        if 'ASP2' not in sample:
            asp2 = asp1
            asp2nma_1 = asp1nma_1
            asp2nms_1 = asp1nms_1
            asp2nma_2 = asp1nma_2
            asp2nms_2 = asp1nms_2
        if 'RSP2' not in sample:
            rsp2 = rsp1
            rsp2nma_1 = rsp1nma_1
            rsp2nms_1 = rsp1nms_1
            rsp2nma_2 = rsp1nma_2
            rsp2nms_2 = rsp1nms_2
        if 'NSP2' not in sample:
            nsp2 = nsp1
            nsp2hq_1 = nsp1hq_1
            nsp2hq_2 = nsp1hq_2
            nsp2nma_1 = nsp1nma_1
            nsp2nma_2 = nsp1nma_2
            nsp2nms_1 = nsp1nms_1
            nsp2nms_2 = nsp1nms_2

        features['ASP1_OVER_RSP1'], features['ASP2_OVER_RSP2'] = asp1/max(1, asp1+rsp1), asp2/max(1, asp2+rsp2)

        features['ASP1_RSP1_1_NM_Z_SCORE'] = Features.calculate_z_score(asp1nma_1, asp1nms_1, asp1, rsp1nma_1, rsp1nms_1, rsp1)
        features['ASP1_RSP1_2_NM_Z_SCORE'] = Features.calculate_z_score(asp1nma_2, asp1nms_2, asp1, rsp1nma_2, rsp1nms_2, rsp1)
        features['ASP2_RSP2_1_NM_Z_SCORE'] = Features.calculate_z_score(asp2nma_1, asp2nms_1, asp2, rsp2nma_1, rsp2nms_1, rsp2)
        features['ASP2_RSP2_2_NM_Z_SCORE'] = Features.calculate_z_score(asp2nma_2, asp2nms_2, asp2, rsp2nma_2, rsp2nms_2, rsp2)

        features['ASP1_NSP1_1_NM_Z_SCORE'] = Features.calculate_z_score(asp1nma_1, asp1nms_1, asp1, nsp1nma_1, nsp1nms_1, nsp1)
        features['ASP1_NSP1_2_NM_Z_SCORE'] = Features.calculate_z_score(asp1nma_2, asp1nms_2, asp1, nsp1nma_2, nsp1nms_2, nsp1)
        features['ASP2_NSP2_1_NM_Z_SCORE'] = Features.calculate_z_score(asp2nma_1, asp2nms_1, asp2, nsp2nma_1, nsp2nms_1, nsp2)
        features['ASP2_NSP2_2_NM_Z_SCORE'] = Features.calculate_z_score(asp2nma_2, asp2nms_2, asp2, nsp2nma_2, nsp2nms_2, nsp2)

        ssp1 = Features.get_number_value(sample, 'SSP1', 0)
        ssp1hq_1, ssp1hq_2 = Features.get_number_value(sample, 'SSP1HQ', [0, 0])
        ssp1nma_1, ssp1nma_2 = Features.get_number_value(sample, 'SSP1NMA', [Features.NAN, Features.NAN])
        ssp1nms_1, ssp1nms_2 = Features.get_number_value(sample, 'SSP1NMS', [Features.NAN, Features.NAN])
        features['SSP1'] = Features.piecewise_normalise(ssp1, min_disc_pairs, max_disc_pairs)
        features['SSP1HQ_1'] = Features.piecewise_normalise(ssp1hq_1, min_disc_pairs, max_disc_pairs)
        features['SSP1HQ_2'] = Features.piecewise_normalise(ssp1hq_2, min_disc_pairs, max_disc_pairs)
        features['SSP1HQ_1_RATIO'], features['SSP1HQ_2_RATIO'] = ssp1hq_1/max(1, ssp1), ssp1hq_2/max(1, ssp1)
        features['SSP1mQ_1'], features['SSP1mQ_2'] = Features.get_number_value(sample, 'SSP1mQ', [Features.NAN, Features.NAN])
        features['SSP1MQ_1'], features['SSP1MQ_2'] = Features.get_number_value(sample, 'SSP1MQ', [Features.NAN, Features.NAN])

        ssp2 = Features.get_number_value(sample, 'SSP2', 0)
        ssp2hq_1, ssp2hq_2 = Features.get_number_value(sample, 'SSP2HQ', [0, 0])
        ssp2nma_1, ssp2nma_2 = Features.get_number_value(sample, 'SSP2NMA', [Features.NAN, Features.NAN])
        ssp2nms_1, ssp2nms_2 = Features.get_number_value(sample, 'SSP2NMS', [Features.NAN, Features.NAN])
        features['SSP2'] = Features.piecewise_normalise(ssp2, min_disc_pairs, max_disc_pairs)
        features['SSP2HQ_1'] = Features.piecewise_normalise(ssp2hq_1, min_disc_pairs, max_disc_pairs)
        features['SSP2HQ_2'] = Features.piecewise_normalise(ssp2hq_2, min_disc_pairs, max_disc_pairs)
        features['SSP2HQ_1_RATIO'], features['SSP2HQ_2_RATIO'] = ssp2hq_1/max(1, ssp2), ssp2hq_2/max(1, ssp2)
        features['SSP2mQ_1'], features['SSP2mQ_2'] = Features.get_number_value(sample, 'SSP2mQ', [Features.NAN, Features.NAN])
        features['SSP2MQ_1'], features['SSP2MQ_2'] = Features.get_number_value(sample, 'SSP2MQ', [Features.NAN, Features.NAN])

        if svtype_str in ["DEL", "INS"]:
            ssp1nma_2 = Features.NAN
            ssp2nma_1 = Features.NAN
        elif svtype_str == "DUP":
            ssp1nma_1 = Features.NAN
            ssp2nma_2 = Features.NAN

        features['SSP1_RSP1_1_NM_Z_SCORE'] = Features.calculate_z_score(ssp1nma_1, ssp1nms_1, ssp1, rsp1nma_1, rsp1nms_1, rsp1)
        features['SSP1_RSP1_2_NM_Z_SCORE'] = Features.calculate_z_score(ssp1nma_2, ssp1nms_2, ssp1, rsp1nma_2, rsp1nms_2, rsp1)
        features['SSP2_RSP2_1_NM_Z_SCORE'] = Features.calculate_z_score(ssp2nma_1, ssp2nms_1, ssp2, rsp2nma_1, rsp2nms_1, rsp2)
        features['SSP2_RSP2_2_NM_Z_SCORE'] = Features.calculate_z_score(ssp2nma_2, ssp2nms_2, ssp2, rsp2nma_2, rsp2nms_2, rsp2)

        axr1, axr2 = Features.get_number_value(sample, 'AXR', [0, 0], median_depth*max_is)
        axr1hq, axr2hq = Features.get_number_value(sample, 'AXRHQ', [0, 0], median_depth*max_is)
        features['AXR1'], features['AXR2'] = axr1, axr2
        features['AXR1HQ'], features['AXR2HQ'] = axr1hq, axr2hq

        exl1 = Features.get_number_value(sample, 'EXL', Features.NAN)
        exl2 = Features.get_number_value(sample, 'EXL2', Features.NAN)
        features['MEXL'] = max(exl1/(max_is+read_len), exl2/(max_is+read_len))
        features['mEXL'] = min(exl1/(max_is+read_len), exl2/(max_is+read_len))
        features['EXL'] = features['MEXL'] + features['mEXL']

        exas1 = Features.get_number_value(sample, 'EXAS', 0)
        exas2 = Features.get_number_value(sample, 'EXAS2', 0)
        exrs1 = Features.get_number_value(sample, 'EXRS', 0)
        exrs2 = Features.get_number_value(sample, 'EXRS2', 0)

        features['EXAS_EXRS_DIFF_TO_LEN'] = (exas1-exrs1+exas2-exrs2)/max(1, edit_distance)

        exss1_1, exss1_2 = Features.get_number_value(sample, 'EXSS', [0, 0])
        features['EXSS1_1'] = exss1_1/(max_is+read_len)
        features['EXSS1_2'] = exss1_2/(max_is+read_len)
        features['EXSS1_RATIO1'] = exss1_1/max(1, exl1)
        features['EXSS1_RATIO2'] = exss1_2/max(1, exl1)

        exss2_1, exss2_2 = Features.get_number_value(sample, 'EXSS2', [0, 0])
        features['EXSS2_1'] = exss2_1/(max_is+read_len)
        features['EXSS2_2'] = exss2_2/(max_is+read_len)
        features['EXSS2_RATIO1'] = exss2_1/max(1, exl2)
        features['EXSS2_RATIO2'] = exss2_2/max(1, exl2)

        exssc1_1, exssc1_2 = Features.get_number_value(sample, 'EXSSC', [Features.NAN, Features.NAN])
        exssc2_1, exssc2_2 = Features.get_number_value(sample, 'EXSSC2', [Features.NAN, Features.NAN])
        exsscia1_1, exsscia1_2 = Features.get_number_value(sample, 'EXSSCIA', [Features.NAN, Features.NAN])
        exsscia2_1, exsscia2_2 = Features.get_number_value(sample, 'EXSSC2IA', [Features.NAN, Features.NAN])
        features['EXSSC1_IA_RATIO'] = (exssc1_1+exssc1_2)/max(1, exsscia1_1+exsscia1_2)
        features['EXSSC2_IA_RATIO'] = (exssc2_1+exssc2_2)/max(1, exsscia2_1+exsscia2_2)
        features['EXSSC1_IA_DIFF'] = (exsscia1_1+exsscia1_2-exssc1_1-exssc1_2)/max(1, exss1_1+exss1_2)
        features['EXSSC2_IA_DIFF'] = (exsscia2_1+exsscia2_2-exssc2_1-exssc2_2)/max(1, exss2_1+exss2_2)

        feature_values = []
        for feature_name in feature_names:
            if feature_name not in features:
                raise RuntimeError(
                    f"Feature '{feature_name}' required for model {model_name} is not produced by features.py."
                )
            feature_values.append(features[feature_name])
        return feature_values

def select_gt(gt1, gt2):
    if gt1 == "./." and gt2 != "./.":
        return gt2
    if gt1 != "./." and gt2 == "./.":
        return gt1
    elif gt1 == gt2:
        return gt1
    else:
        return "./."

def gt_alleles(gt):
    return gt.replace("|", "/").split("/")

def gt_has_alt(gt):
    return "1" in gt_alleles(gt)

def gt_is_known_positive(gt):
    alleles = gt_alleles(gt)
    return "1" in alleles and "." not in alleles

def gt_is_hom_alt(gt):
    alleles = gt_alleles(gt)
    return len(alleles) > 0 and all(allele == "1" for allele in alleles)

def gt_stage_label(gt):
    if gt == "0/1":
        return 0
    if gt == "1/2":
        return 1
    if gt == "1/1":
        return 2
    raise RuntimeError(f"Unsupported GT-stage label: {gt}")

def gt_has_alt_array(gts):
    return np.array([gt_has_alt(gt) for gt in gts])

def gt_is_known_positive_array(gts):
    return np.array([gt_is_known_positive(gt) for gt in gts])

def gt_is_hom_alt_array(gts):
    return np.array([gt_is_hom_alt(gt) for gt in gts])

def read_gt_labels(file_path):
    if not os.path.exists(file_path):
        raise RuntimeError(f"Genotype labels file not found: {file_path}")
    gt_labels = dict()
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 3:
                raise RuntimeError(f"Malformed genotype labels line in {file_path}: {' '.join(fields)}")
            id, gt, exact = fields[0], fields[1], int(fields[2])
            if id not in gt_labels:
                gt_labels[id] = (gt, exact)
            else:
                prev_gt, prev_exact = gt_labels[id]
                gt_labels[id] = (
                    select_gt(prev_gt, gt),
                    max(prev_exact, exact),
                )
    return gt_labels

def load_stats(stats_fname):
    stats = defaultdict(dict)
    with open(stats_fname, 'r') as stats_reader:
        for line in stats_reader:
            sl = line.strip().split()
            stats[sl[0]][sl[1]] = int(sl[2])
    return stats

def get_stat(stats, stat_name, chrom):
    if chrom in stats[stat_name]:
        return stats[stat_name][chrom]
    return stats[stat_name]['.']

# Function to parse the VCF file and extract relevant features using pysam
def parse_vcf(vcf_fname, stats_fname, fp_fname, ignore_gts = False, feature_names_by_model = None, restrict_to_model_name = None, gt_labels = None):
    if not ignore_gts and gt_labels is None:
        gt_labels = read_gt_labels(fp_fname)
    vcf_reader = pysam.VariantFile(vcf_fname)
    stats = load_stats(stats_fname)

    features_by_source, gts_by_source, variant_ids_by_source, exacts_by_source = defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
    for record in vcf_reader.fetch():
        if Features.skips_ml_genotyping(record):
            continue

        model_name = Features.get_model_name(record, get_stat(stats, 'max_is', record.chrom), get_stat(stats, 'read_len', record.chrom))
        if restrict_to_model_name not in (None, "ALL") and model_name != restrict_to_model_name:
            continue
        if ignore_gts:
            gt = "NA"
            exact = "NA"
        else:
            label = gt_labels.get(record.id)
            if label is None:
                raise RuntimeError(f"Missing GT label for record {record.id} in {fp_fname}")
            gt, exact = label
        model_feature_names = None if feature_names_by_model is None else feature_names_by_model.get(model_name)
        feature_values = Features.record_to_features(record, stats, model_feature_names)
        features_by_source[model_name].append(feature_values)
        gts_by_source[model_name].append(gt)
        variant_ids_by_source[model_name].append(Features.generate_id(record))
        exacts_by_source[model_name].append(exact)

    for model_name in features_by_source:
        features_by_source[model_name] = np.array(features_by_source[model_name])
        gts_by_source[model_name] = np.array(gts_by_source[model_name])
        variant_ids_by_source[model_name] = np.array(variant_ids_by_source[model_name])
        exacts_by_source[model_name] = np.array(exacts_by_source[model_name])
    
    return features_by_source, gts_by_source, variant_ids_by_source, exacts_by_source
