from __future__ import division
import os, pysam
from collections import defaultdict
from typing import NamedTuple, Tuple
import numpy as np
import math

class AltReadMetrics(NamedTuple):
    ar1: float
    ar2: float
    ar1c: float
    ar2c: float
    ar1chq: float
    ar2chq: float
    ar1e: float
    ar2e: float
    hp_genotyped: bool
    chrom: str
    start: int
    stop: int
    rr1: float
    rr2: float
    rr1c: float
    rr2c: float
    rr1e: float
    rr2e: float
    oar1_vids: Tuple[str, ...]
    oar2_vids: Tuple[str, ...]
    xaas_xars_diff_to_len: float
    has_extension_evidence: bool
    aas_ars_diff_to_len: float = 0.0
    has_assembly_evidence: bool = False

    def ar(self, bp_idx): return self.ar1 if bp_idx == 0 else self.ar2
    def arc(self, bp_idx): return self.ar1c if bp_idx == 0 else self.ar2c
    def archq(self, bp_idx): return self.ar1chq if bp_idx == 0 else self.ar2chq
    def are(self, bp_idx): return self.ar1e if bp_idx == 0 else self.ar2e
    def oar_vids(self, bp_idx): return self.oar1_vids if bp_idx == 0 else self.oar2_vids

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

    reads_features_names = ['AR1', 'AR1_ADJ', 'AR1C', 'AR1C_ADJ', 'AR1C_RATIO', 'AR1CmQ', 'AR1CMQ', 'AR1CHQ', 'AR1C_HQ_RATIO', 
                            'AR2', 'AR2_ADJ', 'AR2C', 'AR2C_ADJ', 'AR2C_RATIO', 'AR2CmQ', 'AR2CMQ', 'AR2CHQ', 'AR2C_HQ_RATIO', 
                            'AR1E', 'AR1EmQ', 'AR1EMQ', 'AR1E_HQ_RATIO', 'AR1E_RATIO',
                            'AR2E', 'AR2EmQ', 'AR2EMQ', 'AR2E_HQ_RATIO', 'AR2E_RATIO',
                            'AR1CE', 'AR1CEmQ', 'AR1CEMQ', 'AR1CE_HQ_RATIO', 'AR1CE_RATIO',
                            'AR2CE', 'AR2CEmQ', 'AR2CEMQ', 'AR2CE_HQ_RATIO', 'AR2CE_RATIO',
                            'AR1HPMODE', 'AR1CHPMODE', 'AR1CHPIQR', 'AR1HPMODE_AR1CHPMODE_DIFF', 'AR1HPMODE_ALTLEN_DIFF', 'AR1CHPMODE_ALTLEN_DIFF',
                            'AR1CHPmQ', 'AR1CHPMQ', 'AR1CHPAQ', 'AR1CHPSQ', 'AR1HP5PMR', 'AR1HP3PMR',
                            'MAXARCD', 'MAXARCED', 'MAXARED',
                            'RR1', 'RR1C', 'RR1CmQ', 'RR1CMQ', 'RR1C_HQ_RATIO',
                            'RR1E', 'RR1EmQ', 'RR1EMQ', 'RR1E_HQ_RATIO', 'RR1E_RATIO',
                            'RR1CE', 'RR1CEmQ', 'RR1CEMQ', 'RR1CE_HQ_RATIO', 'RR1CE_RATIO',
                            'RR1HPMODE', 'RR1CHPMODE', 'RR1CHPIQR', 'RR1HPMODE_RR1CHPMODE_DIFF', 'RR1HPMODE_REFLEN_DIFF', 'RR1CHPMODE_REFLEN_DIFF',
                            'RR1CHPmQ', 'RR1CHPMQ', 'RR1CHPAQ', 'RR1CHPSQ', 'RR1HP5PMR', 'RR1HP3PMR',
                            'RR2', 'RR2C', 'RR2CmQ', 'RR2CMQ', 'RR2C_HQ_RATIO',
                            'RR2E', 'RR2EmQ', 'RR2EMQ', 'RR2E_HQ_RATIO', 'RR2E_RATIO',
                            'RR2CE', 'RR2CEmQ', 'RR2CEMQ', 'RR2CE_HQ_RATIO', 'RR2CE_RATIO',
                            'MAXRRCD', 'MAXRRCED', 'MAXRRED',
                            'OAR1', 'OAR2', 'OAR1ALL', 'OAR2ALL', 'OAR1_ALL_RATIO', 'OAR2_ALL_RATIO',
                            'OTHER_HP_GENOTYPED',
                            'OTHER1_MIN_ARCHQ', 'OTHER2_MIN_ARCHQ',
                            'OTHER1_MIN_AR_OVER_NAR', 'OTHER2_MIN_AR_OVER_NAR',
                            'OTHER1_MIN_ARC_OVER_NARC', 'OTHER2_MIN_ARC_OVER_NARC',
                            'OTHER1_MIN_ARE_OVER_NARE', 'OTHER2_MIN_ARE_OVER_NARE',
                            'OTHER1_AAS_ARS_DIFF_TO_LEN', 'OTHER2_AAS_ARS_DIFF_TO_LEN',
                            'OTHER1_XAAS_XARS_DIFF_TO_LEN', 'OTHER2_XAAS_XARS_DIFF_TO_LEN',
                            'OAR1C', 'OAR2C', 'OAR1CHQ', 'OAR2CHQ', 'OAR1C_HQ_RATIO', 'OAR2C_HQ_RATIO', 'OAR1E', 'OAR2E',
                            'NAR1', 'NAR2', 'NAR1C', 'NAR2C', 'NAR1CHQ', 'NAR2CHQ', 'NAR1C_HQ_RATIO', 'NAR2C_HQ_RATIO',
                            'NAR1CE', 'NAR2CE', 'NAR1E', 'NAR2E',
                            'ER1_DEVIATION', 'ER2_DEVIATION', 'ER_HQ_RATIO',
                            'AR1CMSPAN_1', 'AR1CMSPAN_2', 'AR1CMHQSPAN_1', 'AR1CMHQSPAN_2',
                            'AR2CMSPAN_1', 'AR2CMSPAN_2', 'AR2CMHQSPAN_1', 'AR2CMHQSPAN_2',
                            'RR1CMSPAN_1', 'RR1CMSPAN_2', 'RR1CMHQSPAN_1', 'RR1CMHQSPAN_2',
                            'RR2CMSPAN_1', 'RR2CMSPAN_2', 'RR2CMHQSPAN_1', 'RR2CMHQSPAN_2',
                            'AR1_RR1_CAS_Z_SCORE', 'AR2_RR2_CAS_Z_SCORE', 
                            'AR1_OVER_RR1', 'AR2_OVER_RR2', 'AR1C_OVER_RR1C', 'AR2C_OVER_RR2C',
                            'AR1CE_OVER_RR1CE', 'AR2CE_OVER_RR2CE', 'AR1E_OVER_RR1E', 'AR2E_OVER_RR2E',
                            'AR1_OVER_OAR1', 'AR2_OVER_OAR2', 'AR1C_OVER_OAR1C', 'AR2C_OVER_OAR2C',
                            'AR1CE_OVER_OAR1E', 'AR2CE_OVER_OAR2E', 'AR1E_OVER_OAR1E', 'AR2E_OVER_OAR2E',
                            'OAR1_OVER_NAR1', 'OAR2_OVER_NAR2', 'OAR1C_OVER_NAR1C', 'OAR2C_OVER_NAR2C',
                            'OAR1E_OVER_NAR1CE', 'OAR2E_OVER_NAR2CE', 'OAR1E_OVER_NAR1E', 'OAR2E_OVER_NAR2E',
                            'OAR1_OVER_TOTAL1', 'OAR2_OVER_TOTAL2', 'OAR1C_OVER_TOTAL1C', 'OAR2C_OVER_TOTAL2C',
                            'OAR1E_OVER_TOTAL1CE', 'OAR2E_OVER_TOTAL2CE', 'OAR1E_OVER_TOTAL1E', 'OAR2E_OVER_TOTAL2E',
                            'ORR1_RATIO', 'ORR2_RATIO', 'ORR1C_RATIO', 'ORR2C_RATIO', 'ORR1E_RATIO', 'ORR2E_RATIO',
                            'AR1_OVER_NAR1', 'AR2_OVER_NAR2', 'AR1C_OVER_NAR1C', 'AR2C_OVER_NAR2C',
                            'AR1CE_OVER_NAR1CE', 'AR2CE_OVER_NAR2CE', 'AR1E_OVER_NAR1E', 'AR2E_OVER_NAR2E']

    fmt_features_names = [  'ASS1_1', 'ASS1_2', 'ASS2_1', 'ASS2_2',
                            'ASS1_RATIO1', 'ASS1_RATIO2', 'ASS2_RATIO1', 'ASS2_RATIO2',
                            'AAS_ARS_DIFF_TO_LEN',
                            'ASSC1_IA_RATIO', 'ASSC2_IA_RATIO', 'ASSC1_IA_DIFF', 'ASSC2_IA_DIFF',
                            'MAL', 'mAL', 'AL',
                            'AXR1', 'AXR2', 'AXR1HQ', 'AXR2HQ',
                            'XASS1_1', 'XASS1_2', 'XASS2_1', 'XASS2_2',
                            'XASS1_RATIO1', 'XASS1_RATIO2', 'XASS2_RATIO1', 'XASS2_RATIO2',
                            'XAAS_XARS_DIFF_TO_LEN',
                            'XASSC1_IA_RATIO', 'XASSC2_IA_RATIO', 'XASSC1_IA_DIFF', 'XASSC2_IA_DIFF',
                            'MXAL', 'mXAL', 'XAL',
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
            if 'XAL' not in sample:
                svtype_str += "_NOEXL"
        elif svtype_str == "DUP" and Features.get_svlen(record) > read_len-30:
            svtype_str += "_LARGE"
            if 'XAL' not in sample:
                svtype_str += "_NOEXL"

        return svtype_str

    def get_number_value(info, key, default, norm_factor = 1.0):
        if key in info:
            v = info[key]
        else:
            v = default
        if isinstance(v, list) or isinstance(v, tuple):
            if isinstance(default, (list, tuple)):
                return [float(default[i] if x is None else x)/norm_factor for i, x in enumerate(v)]
            return [float(default if x is None else x)/norm_factor for x in v]
        else:
            if v is None:
                v = default
            return float(v)/norm_factor

    def exact_read_ratio(exact_reads, reads):
        return exact_reads/max(1, reads)

    def consistent_exact_read_ratio(consistent_exact_reads, consistent_reads):
        return consistent_exact_reads/max(1, consistent_reads)

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

    def normalize_sv_id(sv_id):
        sv_id = str(sv_id)
        return sv_id[:-4] if sv_id.endswith('_DUP') else sv_id

    def get_oar_vids(sample, key):
        vids = []
        for value in Features.get_string_list_value(sample, key, []):
            for vid in value.split(','):
                if vid and vid != '.':
                    vids.append(Features.normalize_sv_id(vid))
        return tuple(vids)

    def select_alt_read_metrics_for_oar_vids(sample, key, bp_idx, alt_reads_by_vid):
        if alt_reads_by_vid is None:
            return 0, 0, 0, 0, 0, 0, 0, 0, 0, Features.NAN, Features.NAN

        candidates = []
        for vid in Features.get_oar_vids(sample, key):
            for values in alt_reads_by_vid.get(vid, []):
                candidates.append((values, vid))

        if not candidates:
            return 0, 0, 0, 0, 0, 0, 0, 0, 0, Features.NAN, Features.NAN

        # Select one source variant by AR, then ARC, then ARE. ARCHQ is carried from that variant but does not affect selection.
        selected, _ = min(candidates, key=lambda x: (-x[0].ar(bp_idx), -x[0].arc(bp_idx), -x[0].are(bp_idx), x[1]))
        selected_oars = []
        for oar_bp_idx in range(2):
            oar_vids = selected.oar_vids(oar_bp_idx)
            oar_candidates = [(values, vid) for vid in oar_vids for values in alt_reads_by_vid.get(vid, [])]
            selected_oars.append(min(oar_candidates, key=lambda x: (-x[0].ar(oar_bp_idx), -x[0].arc(oar_bp_idx), -x[0].are(oar_bp_idx), x[1]))[0] if oar_candidates else None)
        oar1, oar2 = selected_oars
        min_ar_over_nar = min(selected.ar1/max(1, selected.ar1 + (oar1.ar1 if oar1 else 0) + selected.rr1), selected.ar2/max(1, selected.ar2 + (oar2.ar2 if oar2 else 0) + selected.rr2))
        min_arc_over_narc = min(selected.ar1c/max(1, selected.ar1c + (oar1.ar1c if oar1 else 0) + selected.rr1c), selected.ar2c/max(1, selected.ar2c + (oar2.ar2c if oar2 else 0) + selected.rr2c))
        min_are_over_nare = min(selected.ar1e/max(1, selected.ar1e + (oar1.ar1e if oar1 else 0) + selected.rr1e), selected.ar2e/max(1, selected.ar2e + (oar2.ar2e if oar2 else 0) + selected.rr2e))
        xaas_xars_diff_to_len = selected.xaas_xars_diff_to_len if selected.has_extension_evidence else Features.NAN
        aas_ars_diff_to_len = selected.aas_ars_diff_to_len if selected.has_assembly_evidence else Features.NAN
        return selected.ar(bp_idx), selected.arc(bp_idx), selected.archq(bp_idx), selected.are(bp_idx), int(selected.hp_genotyped), min(selected.ar1chq, selected.ar2chq), min_ar_over_nar, min_arc_over_narc, min_are_over_nare, xaas_xars_diff_to_len, aas_ars_diff_to_len

    def consensus_alt_ref_score_diff_to_len(record, prefix):
        sample = record.samples[0]
        aas1 = Features.get_number_value(sample, prefix+'AAS', 0)
        aas2 = Features.get_number_value(sample, prefix+'AAS2', 0)
        ars1 = Features.get_number_value(sample, prefix+'ARS', 0)
        ars2 = Features.get_number_value(sample, prefix+'ARS2', 0)
        svinslen = len(Features.get_svinsseq(record))
        edit_distance = Features.get_edit_distance(record, svinslen)
        return (aas1-ars1+aas2-ars2)/max(1, edit_distance)

    def aas_ars_diff_to_len(record):
        return Features.consensus_alt_ref_score_diff_to_len(record, '')

    def xaas_xars_diff_to_len(record):
        return Features.consensus_alt_ref_score_diff_to_len(record, 'X')

    def has_assembly_evidence(record):
        sample = record.samples[0]
        return 'AL' in sample or 'AL2' in sample

    def has_extension_evidence(record):
        sample = record.samples[0]
        return 'XAL' in sample or 'XAL2' in sample

    def add_consensus_alignment_features(features, sample, prefix, max_is, read_len, edit_distance):
        al1 = Features.get_number_value(sample, prefix+'AL', Features.NAN)
        al2 = Features.get_number_value(sample, prefix+'AL2', Features.NAN)
        length_normalisation_factor = max_is+read_len if prefix == 'X' else read_len
        al_normalisation_factors = [length_normalisation_factor, length_normalisation_factor]
        if prefix == '':
            al_normalisation_factors = Features.get_number_value(sample, 'MFAL', [read_len, read_len])
        normalised_als = [al/al_normalisation_factors[i] for i, al in enumerate((al1, al2)) if math.isfinite(al)]
        if normalised_als:
            features['M'+prefix+'AL'] = max(normalised_als)
            features['m'+prefix+'AL'] = min(normalised_als) if len(normalised_als) == 2 else Features.NAN
            features[prefix+'AL'] = sum(normalised_als)
        else:
            features['M'+prefix+'AL'] = Features.NAN
            features['m'+prefix+'AL'] = Features.NAN
            features[prefix+'AL'] = Features.NAN

        aas1 = Features.get_number_value(sample, prefix+'AAS', 0)
        aas2 = Features.get_number_value(sample, prefix+'AAS2', 0)
        ars1 = Features.get_number_value(sample, prefix+'ARS', 0)
        ars2 = Features.get_number_value(sample, prefix+'ARS2', 0)
        features[prefix+'AAS_'+prefix+'ARS_DIFF_TO_LEN'] = (aas1-ars1+aas2-ars2)/max(1, edit_distance)

        ass1_1, ass1_2 = Features.get_number_value(sample, prefix+'ASS', [Features.NAN, Features.NAN])
        features[prefix+'ASS1_1'] = ass1_1/al_normalisation_factors[0]
        features[prefix+'ASS1_2'] = ass1_2/al_normalisation_factors[0]
        features[prefix+'ASS1_RATIO1'] = ass1_1/max(1, al1)
        features[prefix+'ASS1_RATIO2'] = ass1_2/max(1, al1)

        ass2_1, ass2_2 = Features.get_number_value(sample, prefix+'ASS2', [Features.NAN, Features.NAN])
        features[prefix+'ASS2_1'] = ass2_1/al_normalisation_factors[1]
        features[prefix+'ASS2_2'] = ass2_2/al_normalisation_factors[1]
        features[prefix+'ASS2_RATIO1'] = ass2_1/max(1, al2)
        features[prefix+'ASS2_RATIO2'] = ass2_2/max(1, al2)

        assc1_1, assc1_2 = Features.get_number_value(sample, prefix+'ASSC', [Features.NAN, Features.NAN])
        assc2_1, assc2_2 = Features.get_number_value(sample, prefix+'ASSC2', [Features.NAN, Features.NAN])
        asscia1_1, asscia1_2 = Features.get_number_value(sample, prefix+'ASSCIA', [Features.NAN, Features.NAN])
        asscia2_1, asscia2_2 = Features.get_number_value(sample, prefix+'ASSC2IA', [Features.NAN, Features.NAN])
        features[prefix+'ASSC1_IA_RATIO'] = (assc1_1+assc1_2)/max(1, asscia1_1+asscia1_2)
        features[prefix+'ASSC2_IA_RATIO'] = (assc2_1+assc2_2)/max(1, asscia2_1+asscia2_2)
        features[prefix+'ASSC1_IA_DIFF'] = (asscia1_1+asscia1_2-assc1_1-assc1_2)/max(1, ass1_1+ass1_2)
        features[prefix+'ASSC2_IA_DIFF'] = (asscia2_1+asscia2_2-assc2_1-assc2_2)/max(1, ass2_1+ass2_2)

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

    def record_to_features(record, stats, feature_names = None, alt_reads_by_vid = None):
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
            features['INS_SEQ_COV_SUFFIX_LEN'] = (len(svinsseq)-i-1)/len(svinsseq)

        if 'EARF' in sample:
            exp_alt_reads_freq1, exp_alt_reads_freq2 = sample['EARF']
        if exp_alt_reads_freq1 is None:
            exp_alt_reads_freq1 = Features.NAN
        if exp_alt_reads_freq2 is None:
            exp_alt_reads_freq2 = Features.NAN
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
        ar1ce = Features.get_number_value(sample, 'AR1CE', 0)
        ar1cehq = Features.get_number_value(sample, 'AR1CEHQ', 0)
        ar1ce_min_mq = Features.get_number_value(sample, 'AR1CEmQ', Features.NAN)
        ar1ce_max_mq = Features.get_number_value(sample, 'AR1CEMQ', Features.NAN)
        ar1e = Features.get_number_value(sample, 'AR1E', 0)
        ar1ehq = Features.get_number_value(sample, 'AR1EHQ', 0)
        ar1e_min_mq = Features.get_number_value(sample, 'AR1EmQ', Features.NAN)
        ar1e_max_mq = Features.get_number_value(sample, 'AR1EMQ', Features.NAN)
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
        features['AR1C_RATIO'] = ar1c/max(1, ar1)
        features['AR1CmQ'] = Features.get_number_value(sample, 'AR1CmQ', Features.NAN)
        features['AR1CMQ'] = Features.get_number_value(sample, 'AR1CMQ', Features.NAN)
        features['AR1CHQ'] = Features.piecewise_normalise(arc1hq, min_depth, max_depth)
        features['AR1C_HQ_RATIO'] = arc1hq/max(1, ar1c)
        features['AR1E'] = Features.piecewise_normalise(ar1e, min_depth, max_depth)
        features['AR1EmQ'] = ar1e_min_mq
        features['AR1EMQ'] = ar1e_max_mq
        features['AR1E_HQ_RATIO'] = ar1ehq/max(1, ar1e)
        features['AR1E_RATIO'] = Features.exact_read_ratio(ar1e, ar1)
        features['AR1CE'] = Features.piecewise_normalise(ar1ce, min_depth, max_depth)
        features['AR1CEmQ'] = ar1ce_min_mq
        features['AR1CEMQ'] = ar1ce_max_mq
        features['AR1CE_HQ_RATIO'] = ar1cehq/max(1, ar1ce)
        features['AR1CE_RATIO'] = Features.consistent_exact_read_ratio(ar1ce, ar1c)
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
        arc2hq = Features.get_number_value(sample, 'AR2CHQ', 0)
        ar2e = Features.get_number_value(sample, 'AR2E', 0)
        ar2ehq = Features.get_number_value(sample, 'AR2EHQ', 0)
        ar2e_min_mq = Features.get_number_value(sample, 'AR2EmQ', Features.NAN)
        ar2e_max_mq = Features.get_number_value(sample, 'AR2EMQ', Features.NAN)
        ar2ce = Features.get_number_value(sample, 'AR2CE', 0)
        ar2cehq = Features.get_number_value(sample, 'AR2CEHQ', 0)
        ar2ce_min_mq = Features.get_number_value(sample, 'AR2CEmQ', Features.NAN)
        ar2ce_max_mq = Features.get_number_value(sample, 'AR2CEMQ', Features.NAN)
        ar2cas = Features.get_number_value(sample, 'AR2CAS', Features.NAN)
        ar2css = Features.get_number_value(sample, 'AR2CSS', Features.NAN)
        has_ar2 = 'AR2' in sample
        if not has_ar2:
            ar2 = ar1
            ar2c = ar1c
            arc2hq = arc1hq
            ar2ce = ar1ce
            ar2cehq = ar1cehq
            ar2ce_min_mq = ar1ce_min_mq
            ar2ce_max_mq = ar1ce_max_mq
            ar2e = ar1e
            ar2ehq = ar1ehq
            ar2e_min_mq = ar1e_min_mq
            ar2e_max_mq = ar1e_max_mq
            ar2cas = ar1cas
            ar2css = ar1css
            ar2_adj = ar1_adj
            ar2c_adj = ar1c_adj
        else:
            ar2_adj = ar2
            ar2c_adj = ar2c
            if exp_alt_reads_freq2 > 0:
                ar2_adj = ar2/exp_alt_reads_freq2
                ar2c_adj = ar2c/exp_alt_reads_freq2
        features['AR2'] = Features.piecewise_normalise(ar2, min_depth, max_depth)
        features['AR2C'] = Features.piecewise_normalise(ar2c, min_depth, max_depth)
        features['AR2_ADJ'] = Features.piecewise_normalise(ar2_adj, min_depth, max_depth)
        features['AR2C_ADJ'] = Features.piecewise_normalise(ar2c_adj, min_depth, max_depth)
        features['AR2C_RATIO'] = ar2c/max(1, ar2)
        features['AR2CmQ'] = Features.get_number_value(sample, 'AR2CmQ', Features.NAN)
        features['AR2CMQ'] = Features.get_number_value(sample, 'AR2CMQ', Features.NAN)
        features['AR2CHQ'] = Features.piecewise_normalise(arc2hq, min_depth, max_depth)
        features['AR2C_HQ_RATIO'] = arc2hq/max(1, ar2c)
        features['AR2E'] = Features.piecewise_normalise(ar2e, min_depth, max_depth)
        features['AR2EmQ'] = ar2e_min_mq
        features['AR2EMQ'] = ar2e_max_mq
        features['AR2E_HQ_RATIO'] = ar2ehq/max(1, ar2e)
        features['AR2E_RATIO'] = Features.exact_read_ratio(ar2e, ar2)
        features['AR2CE'] = Features.piecewise_normalise(ar2ce, min_depth, max_depth)
        features['AR2CEmQ'] = ar2ce_min_mq
        features['AR2CEMQ'] = ar2ce_max_mq
        features['AR2CE_HQ_RATIO'] = ar2cehq/max(1, ar2ce)
        features['AR2CE_RATIO'] = Features.consistent_exact_read_ratio(ar2ce, ar2c)
        features['AR2CMSPAN_1'], features['AR2CMSPAN_2'] = Features.get_number_value(sample, 'AR2CMSPAN', [0, 0], max_is)
        features['AR2CMHQSPAN_1'], features['AR2CMHQSPAN_2'] = Features.get_number_value(sample, 'AR2CMHQSPAN', [0, 0], max_is)
        if not has_ar2:
            features['AR2CmQ'], features['AR2CMQ'] = features['AR1CmQ'], features['AR1CMQ']
            features['AR2CMSPAN_1'], features['AR2CMSPAN_2'] = features['AR1CMSPAN_1'], features['AR1CMSPAN_2']
            features['AR2CMHQSPAN_1'], features['AR2CMHQSPAN_2'] = features['AR1CMHQSPAN_1'], features['AR1CMHQSPAN_2']

        ar1cf = Features.get_number_value(sample, 'AR1CF', 0, max(1, ar1c))
        ar1cr = Features.get_number_value(sample, 'AR1CR', 0, max(1, ar1c))
        ar2cf = Features.get_number_value(sample, 'AR2CF', 0, max(1, ar2c))
        ar2cr = Features.get_number_value(sample, 'AR2CR', 0, max(1, ar2c))
        if not has_ar2: ar2cf, ar2cr = ar1cf, ar1cr
        features['ARCF'] = ar1cf + ar2cf
        features['ARCR'] = ar1cr + ar2cr
        features['MAXARCD'] = max(features['ARCF'], features['ARCR'])

        ar1ef = Features.get_number_value(sample, 'AR1CEF', 0, max(1, ar1ce))
        ar1er = Features.get_number_value(sample, 'AR1CER', 0, max(1, ar1ce))
        ar2ef = Features.get_number_value(sample, 'AR2CEF', 0, max(1, ar2ce))
        ar2er = Features.get_number_value(sample, 'AR2CER', 0, max(1, ar2ce))
        if not has_ar2: ar2ef, ar2er = ar1ef, ar1er
        features['ARCEF'] = ar1ef + ar2ef
        features['ARCER'] = ar1er + ar2er
        features['MAXARCED'] = max(features['ARCEF'], features['ARCER'])

        ar1_exact_fwd = Features.get_number_value(sample, 'AR1EF', 0, max(1, ar1e))
        ar1_exact_rev = Features.get_number_value(sample, 'AR1ER', 0, max(1, ar1e))
        ar2_exact_fwd = Features.get_number_value(sample, 'AR2EF', 0, max(1, ar2e))
        ar2_exact_rev = Features.get_number_value(sample, 'AR2ER', 0, max(1, ar2e))
        if not has_ar2: ar2_exact_fwd, ar2_exact_rev = ar1_exact_fwd, ar1_exact_rev
        features['AREF'] = ar1_exact_fwd + ar2_exact_fwd
        features['ARER'] = ar1_exact_rev + ar2_exact_rev
        features['MAXARED'] = max(features['AREF'], features['ARER'])

        oar1, oar1c, oar1chq, oar1e, oar1_hp_genotyped, other1_min_archq, other1_min_ar_over_nar, other1_min_arc_over_narc, other1_min_are_over_nare, other1_xaas_xars_diff_to_len, other1_aas_ars_diff_to_len = Features.select_alt_read_metrics_for_oar_vids(sample, 'OAR1VID', 0, alt_reads_by_vid)
        oar1all = Features.get_number_value(sample, 'OAR1ALL', 0)
        oar2, oar2c, oar2chq, oar2e, oar2_hp_genotyped, other2_min_archq, other2_min_ar_over_nar, other2_min_arc_over_narc, other2_min_are_over_nare, other2_xaas_xars_diff_to_len, other2_aas_ars_diff_to_len = Features.select_alt_read_metrics_for_oar_vids(sample, 'OAR2VID', 1, alt_reads_by_vid)
        oar2all = Features.get_number_value(sample, 'OAR2ALL', 0)
        if not has_ar2:
            oar2 = oar1
            oar2all = oar1all
            oar2c = oar1c
            oar2chq = oar1chq
            oar2e = oar1e
            oar2_hp_genotyped = oar1_hp_genotyped
            other2_min_archq = other1_min_archq
            other2_min_ar_over_nar = other1_min_ar_over_nar
            other2_min_arc_over_narc = other1_min_arc_over_narc
            other2_min_are_over_nare = other1_min_are_over_nare
            other2_xaas_xars_diff_to_len = other1_xaas_xars_diff_to_len
            other2_aas_ars_diff_to_len = other1_aas_ars_diff_to_len

        features['OAR1'] = Features.piecewise_normalise(oar1, min_depth, max_depth)
        features['OAR1ALL'] = Features.piecewise_normalise(oar1all, min_depth, max_depth)
        features['OAR1_ALL_RATIO'] = oar1/max(1, oar1all)
        features['OTHER_HP_GENOTYPED'] = max(oar1_hp_genotyped, oar2_hp_genotyped)
        features['OTHER1_MIN_ARC_OVER_NARC'] = other1_min_arc_over_narc
        features['OTHER1_MIN_ARCHQ'] = Features.piecewise_normalise(other1_min_archq, min_depth, max_depth)
        features['OTHER1_MIN_AR_OVER_NAR'] = other1_min_ar_over_nar
        features['OTHER1_MIN_ARE_OVER_NARE'] = other1_min_are_over_nare
        features['OTHER1_AAS_ARS_DIFF_TO_LEN'] = other1_aas_ars_diff_to_len
        features['OTHER1_XAAS_XARS_DIFF_TO_LEN'] = other1_xaas_xars_diff_to_len
        features['OAR1C'] = Features.piecewise_normalise(oar1c, min_depth, max_depth)
        features['OAR1CHQ'] = Features.piecewise_normalise(oar1chq, min_depth, max_depth)
        features['OAR1C_HQ_RATIO'] = oar1chq/max(1, oar1c)
        features['OAR1E'] = Features.piecewise_normalise(oar1e, min_depth, max_depth)
        features['OAR2'] = Features.piecewise_normalise(oar2, min_depth, max_depth)
        features['OAR2ALL'] = Features.piecewise_normalise(oar2all, min_depth, max_depth)
        features['OAR2_ALL_RATIO'] = oar2/max(1, oar2all)
        features['OTHER2_MIN_ARC_OVER_NARC'] = other2_min_arc_over_narc
        features['OTHER2_MIN_ARCHQ'] = Features.piecewise_normalise(other2_min_archq, min_depth, max_depth)
        features['OTHER2_MIN_AR_OVER_NAR'] = other2_min_ar_over_nar
        features['OTHER2_MIN_ARE_OVER_NARE'] = other2_min_are_over_nare
        features['OTHER2_AAS_ARS_DIFF_TO_LEN'] = other2_aas_ars_diff_to_len
        features['OTHER2_XAAS_XARS_DIFF_TO_LEN'] = other2_xaas_xars_diff_to_len
        features['OAR2C'] = Features.piecewise_normalise(oar2c, min_depth, max_depth)
        features['OAR2CHQ'] = Features.piecewise_normalise(oar2chq, min_depth, max_depth)
        features['OAR2C_HQ_RATIO'] = oar2chq/max(1, oar2c)
        features['OAR2E'] = Features.piecewise_normalise(oar2e, min_depth, max_depth)

        orr1 = Features.get_number_value(sample, 'ORR1', 0)
        orr1c = Features.get_number_value(sample, 'ORR1C', 0)
        orr1e = Features.get_number_value(sample, 'ORR1E', 0)
        orr2 = Features.get_number_value(sample, 'ORR2', 0)
        orr2c = Features.get_number_value(sample, 'ORR2C', 0)
        orr2e = Features.get_number_value(sample, 'ORR2E', 0)
        if 'ORR2' not in sample:
            orr2 = orr1
            orr2c = orr1c
            orr2e = orr1e

        rr1 = Features.get_number_value(sample, 'RR1', 0)
        rr1c = Features.get_number_value(sample, 'RR1C', 0)
        rr1chq = Features.get_number_value(sample, 'RR1CHQ', 0)
        rr1e = Features.get_number_value(sample, 'RR1E', 0)
        rr1ehq = Features.get_number_value(sample, 'RR1EHQ', 0)
        rr1e_min_mq = Features.get_number_value(sample, 'RR1EmQ', Features.NAN)
        rr1e_max_mq = Features.get_number_value(sample, 'RR1EMQ', Features.NAN)
        rr1ce = Features.get_number_value(sample, 'RR1CE', 0)
        rr1cehq = Features.get_number_value(sample, 'RR1CEHQ', 0)
        rr1ce_min_mq = Features.get_number_value(sample, 'RR1CEmQ', Features.NAN)
        rr1ce_max_mq = Features.get_number_value(sample, 'RR1CEMQ', Features.NAN)
        rr2 = Features.get_number_value(sample, 'RR2', 0)
        rr2c = Features.get_number_value(sample, 'RR2C', 0)
        rr2chq = Features.get_number_value(sample, 'RR2CHQ', 0)
        rr2e = Features.get_number_value(sample, 'RR2E', 0)
        rr2ehq = Features.get_number_value(sample, 'RR2EHQ', 0)
        rr2e_min_mq = Features.get_number_value(sample, 'RR2EmQ', Features.NAN)
        rr2e_max_mq = Features.get_number_value(sample, 'RR2EMQ', Features.NAN)
        rr2ce = Features.get_number_value(sample, 'RR2CE', 0)
        rr2cehq = Features.get_number_value(sample, 'RR2CEHQ', 0)
        rr2ce_min_mq = Features.get_number_value(sample, 'RR2CEmQ', Features.NAN)
        rr2ce_max_mq = Features.get_number_value(sample, 'RR2CEMQ', Features.NAN)

        rr1cas = Features.get_number_value(sample, 'RR1CAS', Features.NAN)
        rr1css = Features.get_number_value(sample, 'RR1CSS', Features.NAN)
        features['RR1'] = Features.piecewise_normalise(rr1, min_depth, max_depth)
        features['RR1C'] = Features.piecewise_normalise(rr1c, min_depth, max_depth)
        features['RR1CmQ'] = Features.get_number_value(sample, 'RR1CmQ', Features.NAN)
        features['RR1CMQ'] = Features.get_number_value(sample, 'RR1CMQ', Features.NAN)
        features['RR1C_HQ_RATIO'] = rr1chq/max(1, rr1c)
        features['RR1E'] = Features.piecewise_normalise(rr1e, min_depth, max_depth)
        features['RR1EmQ'] = rr1e_min_mq
        features['RR1EMQ'] = rr1e_max_mq
        features['RR1E_HQ_RATIO'] = rr1ehq/max(1, rr1e)
        features['RR1E_RATIO'] = Features.exact_read_ratio(rr1e, rr1)
        features['RR1CE'] = Features.piecewise_normalise(rr1ce, min_depth, max_depth)
        features['RR1CEmQ'] = rr1ce_min_mq
        features['RR1CEMQ'] = rr1ce_max_mq
        features['RR1CE_HQ_RATIO'] = rr1cehq/max(1, rr1ce)
        features['RR1CE_RATIO'] = Features.consistent_exact_read_ratio(rr1ce, rr1c)
        features['RR1CMSPAN_1'], features['RR1CMSPAN_2'] = Features.get_number_value(sample, 'RR1CMSPAN', [0, 0], max_is)
        features['RR1CMHQSPAN_1'], features['RR1CMHQSPAN_2'] = Features.get_number_value(sample, 'RR1CMHQSPAN', [0, 0], max_is)

        features['RR1HPMODE'] = Features.get_number_value(sample, 'RR1HPMODE', Features.NAN)
        features['RR1CHPMODE'] = Features.get_number_value(sample, 'RR1CHPMODE', Features.NAN)
        features['RR1HPMODE_RR1CHPMODE_DIFF'] = features['RR1HPMODE'] - features['RR1CHPMODE']
        features['RR1HPMODE_REFLEN_DIFF'] = features['RR1HPMODE'] - features['HP_REF_LEN']
        features['RR1CHPMODE_REFLEN_DIFF'] = features['RR1CHPMODE'] - features['HP_REF_LEN']
        features['RR1CHPIQR'] = Features.get_number_value(sample, 'RR1CHPIQR', Features.NAN)
        features['RR1CHPmQ'] = Features.get_number_value(sample, 'RR1CHPmQ', Features.NAN)
        features['RR1CHPMQ'] = Features.get_number_value(sample, 'RR1CHPMQ', Features.NAN)
        features['RR1CHPAQ'] = Features.get_number_value(sample, 'RR1CHPAQ', Features.NAN)
        features['RR1CHPSQ'] = Features.get_number_value(sample, 'RR1CHPSQ', Features.NAN)
        features['RR1HP5PMR'] = Features.get_number_value(sample, 'RR1HP5PMR', Features.NAN)
        features['RR1HP3PMR'] = Features.get_number_value(sample, 'RR1HP3PMR', Features.NAN)

        rr2cas = Features.get_number_value(sample, 'RR2CAS', Features.NAN)
        rr2css = Features.get_number_value(sample, 'RR2CSS', Features.NAN)
        has_rr2 = 'RR2' in sample
        if not has_rr2:
            rr2 = rr1
            rr2c = rr1c
            rr2chq = rr1chq
            rr2ce = rr1ce
            rr2cehq = rr1cehq
            rr2ce_min_mq = rr1ce_min_mq
            rr2ce_max_mq = rr1ce_max_mq
            rr2e = rr1e
            rr2ehq = rr1ehq
            rr2e_min_mq = rr1e_min_mq
            rr2e_max_mq = rr1e_max_mq
            rr2cas = rr1cas
            rr2css = rr1css
        features['RR2'] = Features.piecewise_normalise(rr2, min_depth, max_depth)
        features['RR2C'] = Features.piecewise_normalise(rr2c, min_depth, max_depth)
        features['RR2CmQ'] = Features.get_number_value(sample, 'RR2CmQ', Features.NAN)
        features['RR2CMQ'] = Features.get_number_value(sample, 'RR2CMQ', Features.NAN)
        features['RR2C_HQ_RATIO'] = rr2chq/max(1, rr2c)
        features['RR2E'] = Features.piecewise_normalise(rr2e, min_depth, max_depth)
        features['RR2EmQ'] = rr2e_min_mq
        features['RR2EMQ'] = rr2e_max_mq
        features['RR2E_HQ_RATIO'] = rr2ehq/max(1, rr2e)
        features['RR2E_RATIO'] = Features.exact_read_ratio(rr2e, rr2)
        features['RR2CE'] = Features.piecewise_normalise(rr2ce, min_depth, max_depth)
        features['RR2CEmQ'] = rr2ce_min_mq
        features['RR2CEMQ'] = rr2ce_max_mq
        features['RR2CE_HQ_RATIO'] = rr2cehq/max(1, rr2ce)
        features['RR2CE_RATIO'] = Features.consistent_exact_read_ratio(rr2ce, rr2c)
        features['RR2CMSPAN_1'], features['RR2CMSPAN_2'] = Features.get_number_value(sample, 'RR2CMSPAN', [0, 0], max_is)
        features['RR2CMHQSPAN_1'], features['RR2CMHQSPAN_2'] = Features.get_number_value(sample, 'RR2CMHQSPAN', [0, 0], max_is)
        if not has_rr2:
            features['RR2CmQ'], features['RR2CMQ'] = features['RR1CmQ'], features['RR1CMQ']
            features['RR2CMSPAN_1'], features['RR2CMSPAN_2'] = features['RR1CMSPAN_1'], features['RR1CMSPAN_2']
            features['RR2CMHQSPAN_1'], features['RR2CMHQSPAN_2'] = features['RR1CMHQSPAN_1'], features['RR1CMHQSPAN_2']

        rr1cf = Features.get_number_value(sample, 'RR1CF', 0)
        rr1cr = Features.get_number_value(sample, 'RR1CR', 0)
        rr2cf = Features.get_number_value(sample, 'RR2CF', 0)
        rr2cr = Features.get_number_value(sample, 'RR2CR', 0)
        if not has_rr2: rr2cf, rr2cr = rr1cf, rr1cr
        rr1cf_ratio = rr1cf/max(1, rr1cf + rr1cr)
        rr1cr_ratio = rr1cr/max(1, rr1cf + rr1cr)
        rr2cf_ratio = rr2cf/max(1, rr2cf + rr2cr)
        rr2cr_ratio = rr2cr/max(1, rr2cf + rr2cr)
        features['RRCF'] = rr1cf_ratio + rr2cf_ratio
        features['RRCR'] = rr1cr_ratio + rr2cr_ratio
        features['MAXRRCD'] = max(features['RRCF'], features['RRCR'])

        rr1ef = Features.get_number_value(sample, 'RR1CEF', 0, max(1, rr1ce))
        rr1er = Features.get_number_value(sample, 'RR1CER', 0, max(1, rr1ce))
        rr2ef = Features.get_number_value(sample, 'RR2CEF', 0, max(1, rr2ce))
        rr2er = Features.get_number_value(sample, 'RR2CER', 0, max(1, rr2ce))
        if not has_rr2: rr2ef, rr2er = rr1ef, rr1er
        features['RRCEF'] = rr1ef + rr2ef
        features['RRCER'] = rr1er + rr2er
        features['MAXRRCED'] = max(features['RRCEF'], features['RRCER'])

        rr1e_fwd = Features.get_number_value(sample, 'RR1EF', 0, max(1, rr1e))
        rr1e_rev = Features.get_number_value(sample, 'RR1ER', 0, max(1, rr1e))
        rr2_exact_fwd = Features.get_number_value(sample, 'RR2EF', 0, max(1, rr2e))
        rr2_exact_rev = Features.get_number_value(sample, 'RR2ER', 0, max(1, rr2e))
        if not has_rr2: rr2_exact_fwd, rr2_exact_rev = rr1e_fwd, rr1e_rev
        features['RREF'] = rr1e_fwd + rr2_exact_fwd
        features['RRER'] = rr1e_rev + rr2_exact_rev
        features['MAXRRED'] = max(features['RREF'], features['RRER'])

        nar1 = rr1 + oar1
        nar1c = rr1c + oar1c
        nar1chq = rr1chq + oar1chq
        nar1e = rr1ce + oar1e
        nar1_exact = rr1e + oar1e
        features['NAR1'] = Features.piecewise_normalise(nar1, min_depth, max_depth)
        features['NAR1C'] = Features.piecewise_normalise(nar1c, min_depth, max_depth)
        features['NAR1CHQ'] = Features.piecewise_normalise(nar1chq, min_depth, max_depth)
        features['NAR1C_HQ_RATIO'] = nar1chq/max(1, nar1c)
        features['NAR1CE'] = Features.piecewise_normalise(nar1e, min_depth, max_depth)
        features['NAR1E'] = Features.piecewise_normalise(nar1_exact, min_depth, max_depth)

        nar2 = rr2 + oar2
        nar2c = rr2c + oar2c
        nar2chq = rr2chq + oar2chq
        nar2e = rr2ce + oar2e
        nar2_exact = rr2e + oar2e
        features['NAR2'] = Features.piecewise_normalise(nar2, min_depth, max_depth)
        features['NAR2C'] = Features.piecewise_normalise(nar2c, min_depth, max_depth)
        features['NAR2CHQ'] = Features.piecewise_normalise(nar2chq, min_depth, max_depth)
        features['NAR2C_HQ_RATIO'] = nar2chq/max(1, nar2c)
        features['NAR2CE'] = Features.piecewise_normalise(nar2e, min_depth, max_depth)
        features['NAR2E'] = Features.piecewise_normalise(nar2_exact, min_depth, max_depth)

        er = Features.get_number_value(sample, 'ER', 0)
        erhq = Features.get_number_value(sample, 'ERHQ', 0)
        er1_total = er+ar1+rr1
        er2_total = er+ar2+rr2
        features['ER1_DEVIATION'] = er/er1_total - (1-exp_alt_reads_freq1) if er1_total > 0 else Features.NAN
        features['ER2_DEVIATION'] = er/er2_total - (1-exp_alt_reads_freq2) if er2_total > 0 else Features.NAN
        features['ER_HQ_RATIO'] = erhq/er if er > 0 else Features.NAN

        features['AR1_RR1_CAS_Z_SCORE'] = Features.calculate_z_score(ar1cas, ar1css, ar1c, rr1cas, rr1css, rr1c)
        features['AR2_RR2_CAS_Z_SCORE'] = Features.calculate_z_score(ar2cas, ar2css, ar2c, rr2cas, rr2css, rr2c)

        features['AR1_OVER_RR1'] = ar1/max(1, ar1+rr1)
        features['AR2_OVER_RR2'] = ar2/max(1, ar2+rr2)
        features['AR1C_OVER_RR1C'] = ar1c/max(1, ar1c+rr1c)
        features['AR2C_OVER_RR2C'] = ar2c/max(1, ar2c+rr2c)
        features['AR1CE_OVER_RR1CE'] = ar1ce/max(1, ar1ce+rr1ce)
        features['AR2CE_OVER_RR2CE'] = ar2ce/max(1, ar2ce+rr2ce)
        features['AR1E_OVER_RR1E'] = ar1e/max(1, ar1e+rr1e)
        features['AR2E_OVER_RR2E'] = ar2e/max(1, ar2e+rr2e)
        features['AR1_OVER_OAR1'] = ar1/max(1, ar1+oar1)
        features['AR2_OVER_OAR2'] = ar2/max(1, ar2+oar2)
        features['AR1C_OVER_OAR1C'] = ar1c/max(1, ar1c+oar1c)
        features['AR2C_OVER_OAR2C'] = ar2c/max(1, ar2c+oar2c)
        features['AR1CE_OVER_OAR1E'] = ar1ce/max(1, ar1ce+oar1e)
        features['AR2CE_OVER_OAR2E'] = ar2ce/max(1, ar2ce+oar2e)
        features['AR1E_OVER_OAR1E'] = ar1e/max(1, ar1e+oar1e)
        features['AR2E_OVER_OAR2E'] = ar2e/max(1, ar2e+oar2e)
        features['OAR1_OVER_NAR1'] = oar1/max(1, nar1)
        features['OAR2_OVER_NAR2'] = oar2/max(1, nar2)
        features['OAR1C_OVER_NAR1C'] = oar1c/max(1, nar1c)
        features['OAR2C_OVER_NAR2C'] = oar2c/max(1, nar2c)
        features['OAR1E_OVER_NAR1CE'] = oar1e/max(1, nar1e)
        features['OAR2E_OVER_NAR2CE'] = oar2e/max(1, nar2e)
        features['OAR1E_OVER_NAR1E'] = oar1e/max(1, nar1_exact)
        features['OAR2E_OVER_NAR2E'] = oar2e/max(1, nar2_exact)
        features['OAR1_OVER_TOTAL1'] = oar1/max(1, nar1+ar1)
        features['OAR2_OVER_TOTAL2'] = oar2/max(1, nar2+ar2)
        features['OAR1C_OVER_TOTAL1C'] = oar1c/max(1, nar1c+ar1c)
        features['OAR2C_OVER_TOTAL2C'] = oar2c/max(1, nar2c+ar2c)
        features['OAR1E_OVER_TOTAL1CE'] = oar1e/max(1, nar1e+ar1ce)
        features['OAR2E_OVER_TOTAL2CE'] = oar2e/max(1, nar2e+ar2ce)
        features['OAR1E_OVER_TOTAL1E'] = oar1e/max(1, nar1_exact+ar1e)
        features['OAR2E_OVER_TOTAL2E'] = oar2e/max(1, nar2_exact+ar2e)
        features['ORR1_RATIO'] = orr1/max(1, rr1)
        features['ORR2_RATIO'] = orr2/max(1, rr2)
        features['ORR1C_RATIO'] = orr1c/max(1, rr1)
        features['ORR2C_RATIO'] = orr2c/max(1, rr2)
        features['ORR1E_RATIO'] = orr1e/max(1, rr1)
        features['ORR2E_RATIO'] = orr2e/max(1, rr2)
        features['AR1_OVER_NAR1'] = ar1/max(1, ar1+nar1)
        features['AR2_OVER_NAR2'] = ar2/max(1, ar2+nar2)
        features['AR1C_OVER_NAR1C'] = ar1c/max(1, ar1c+nar1c)
        features['AR2C_OVER_NAR2C'] = ar2c/max(1, ar2c+nar2c)
        features['AR1CE_OVER_NAR1CE'] = ar1ce/max(1, ar1ce+nar1e)
        features['AR2CE_OVER_NAR2CE'] = ar2ce/max(1, ar2ce+nar2e)
        features['AR1E_OVER_NAR1E'] = ar1e/max(1, ar1e+nar1_exact)
        features['AR2E_OVER_NAR2E'] = ar2e/max(1, ar2e+nar2_exact)

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

        Features.add_consensus_alignment_features(features, sample, '', max_is, read_len, edit_distance)
        Features.add_consensus_alignment_features(features, sample, 'X', max_is, read_len, edit_distance)

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

def gt_is_known_het(gt):
    alleles = gt_alleles(gt)
    return "." not in alleles and "1" in alleles and not gt_is_hom_alt(gt)

def gt_stage_label(gt):
    if gt == "0/1" or gt == "1/2":
        return 0
    if gt == "1/1":
        return 1
    raise RuntimeError(f"Unsupported GT-stage label: {gt}")

def gt_eref_label(gt):
    if gt == "0/1":
        return 1
    if gt == "1/2":
        return 0
    raise RuntimeError(f"Unsupported EREFA-stage label: {gt}")

def gt_has_alt_array(gts):
    return np.array([gt_has_alt(gt) for gt in gts])

def gt_is_known_positive_array(gts):
    return np.array([gt_is_known_positive(gt) for gt in gts])

def gt_is_hom_alt_array(gts):
    return np.array([gt_is_hom_alt(gt) for gt in gts])

def gt_is_known_het_array(gts):
    return np.array([gt_is_known_het(gt) for gt in gts])

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

    alt_reads_by_vid = defaultdict(list)
    for candidate in vcf_reader.fetch():
        candidate_sample = candidate.samples[0]
        ar1 = Features.get_number_value(candidate_sample, 'AR1', 0)
        ar2 = Features.get_number_value(candidate_sample, 'AR2', ar1)
        rr1 = Features.get_number_value(candidate_sample, 'RR1', 0)
        rr2 = Features.get_number_value(candidate_sample, 'RR2', rr1)
        rr1c = Features.get_number_value(candidate_sample, 'RR1C', 0)
        rr2c = Features.get_number_value(candidate_sample, 'RR2C', rr1c)
        rr1e = Features.get_number_value(candidate_sample, 'RR1E', 0)
        rr2e = Features.get_number_value(candidate_sample, 'RR2E', rr1e)
        oar1_vids = Features.get_oar_vids(candidate_sample, 'OAR1VID')
        oar2_vids = Features.get_oar_vids(candidate_sample, 'OAR2VID') if 'AR2' in candidate_sample else oar1_vids
        ar1c = Features.get_number_value(candidate_sample, 'AR1C', 0)
        ar2c = Features.get_number_value(candidate_sample, 'AR2C', ar1c)
        ar1chq = Features.get_number_value(candidate_sample, 'AR1CHQ', 0)
        ar2chq = Features.get_number_value(candidate_sample, 'AR2CHQ', ar1chq)
        ar1e = Features.get_number_value(candidate_sample, 'AR1E', 0)
        ar2e = Features.get_number_value(candidate_sample, 'AR2E', ar1e)
        candidate_xaas_xars_diff_to_len = Features.xaas_xars_diff_to_len(candidate)
        candidate_has_extension_evidence = Features.has_extension_evidence(candidate)
        candidate_aas_ars_diff_to_len = Features.aas_ars_diff_to_len(candidate)
        candidate_has_assembly_evidence = Features.has_assembly_evidence(candidate)
        candidate_id = Features.normalize_sv_id(candidate.id)
        alt_reads_by_vid[candidate_id].append(AltReadMetrics(
            ar1=ar1, ar2=ar2, ar1c=ar1c, ar2c=ar2c, ar1chq=ar1chq, ar2chq=ar2chq, ar1e=ar1e, ar2e=ar2e,
            hp_genotyped=Features.gt_as_homopolymer(candidate), chrom=candidate.chrom, start=candidate.start, stop=candidate.stop,
            rr1=rr1, rr2=rr2, rr1c=rr1c, rr2c=rr2c, rr1e=rr1e, rr2e=rr2e, oar1_vids=oar1_vids, oar2_vids=oar2_vids,
            xaas_xars_diff_to_len=candidate_xaas_xars_diff_to_len, has_extension_evidence=candidate_has_extension_evidence,
            aas_ars_diff_to_len=candidate_aas_ars_diff_to_len, has_assembly_evidence=candidate_has_assembly_evidence))
    vcf_reader.close()
    vcf_reader = pysam.VariantFile(vcf_fname)

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
        feature_values = Features.record_to_features(record, stats, model_feature_names, alt_reads_by_vid)
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
