import os
def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes
def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    allowed template fields - follow python string module:
    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    
    t1w = create_key('sub-{subject}/anat/sub-{subject}_run-01_T1w')
    inv1 = create_key('sub-{subject}/anat/sub-{subject}_run-01_INV1')
    inv2 = create_key('sub-{subject}/anat/sub-{subject}_run-01_INV2')

    func_task01 = create_key('sub-{subject}/func/sub-{subject}_task-semantic_acq-SESB_run-01_bold')
    func_task02 = create_key('sub-{subject}/fmap/sub-{subject}_task-semantic_acq-SESB_dir-PA_run-01_epi')

    func_task03 = create_key('sub-{subject}/func/sub-{subject}_task-semantic_acq-SEMB_run-01_bold')
    func_task04 = create_key('sub-{subject}/fmap/sub-{subject}_task-semantic_acq-SEMB_dir-PA_run-01_epi')

    func_task05 = create_key('sub-{subject}/func/sub-{subject}_task-semantic_acq-MESB_run-01_bold')
    func_task06 = create_key('sub-{subject}/fmap/sub-{subject}_task-semantic_acq-MESB_dir-PA_run-01_epi')

    func_task07 = create_key('sub-{subject}/func/sub-{subject}_task-semantic_acq-MEMB_run-01_bold')
    func_task08 = create_key('sub-{subject}/fmap/sub-{subject}_task-semantic_acq-MEMB_dir-PA_run-01_epi')

    func_task09 = create_key('sub-{subject}/func/sub-{subject}_task-semantic_acq-ptx8ms_run-01_bold')
    func_task10 = create_key('sub-{subject}/fmap/sub-{subject}_task-semantic_acq-ptx8ms_dir-PA_run-01_epi')

    info = {t1w: [], inv1: [], inv2: [], func_task01: [], func_task02: [], func_task03: [], func_task04: [], func_task05: [], func_task06: [], func_task07: [], func_task08: [], func_task09: [], func_task10: []}
    
    for idx, s in enumerate(seqinfo):
       if (s.dim1 == 320) and (s.dim2 == 300) and ('mp2rage_sag_p3_0.75mm' in s.protocol_name):
          if (s.is_derived):
            info[t1w].append(s.series_id)
       if (s.dim1 == 320) and (s.dim2 == 300) and ('mp2rage_sag_p3_0.75mm_INV1' in s.series_description):
            info[inv1].append(s.series_id)
       if (s.dim1 == 320) and (s.dim2 == 300) and ('mp2rage_sag_p3_0.75mm_INV2' in s.series_description):
            info[inv2].append(s.series_id)
       
       if (s.dim4 == 171) and ('Run1_cmrr_mbep2d_bold_25mm_NO_GRAPPA_PF78_SE_MB1' in s.series_description):
            info[func_task01].append(s.series_id)
       if (s.dim4 == 5) and ('Run1_inv_cmrr_mbep2d_bold_25mm_NO_GRAPPA_PF78_SE_MB1' in s.series_description):
            info[func_task02].append(s.series_id)
       if (s.dim4 == 340) and ('Run2_cmrr_mbep2d_bold_25mm_NO_GRAPPA_PF78_SE_MB2' in s.series_description):
            info[func_task03].append(s.series_id)
       if (s.dim4 == 5) and ('Run2_inv_cmrr_mbep2d_bold_25mm_NO_GRAPPA_PF78_SE_MB2' in s.series_description):
            info[func_task04].append(s.series_id)
       if (s.dim4 == 513) and ('Run3_cmrr_mbep2d_bold_25mm_GRAPPA3_PF78_ME3_MB1' in s.series_description):
            info[func_task05].append(s.series_id)
       if (s.dim4 == 15) and ('Run3_inv_cmrr_mbep2d_bold_25mm_GRAPPA3_PF78_ME3_MB1' in s.series_description):
            info[func_task06].append(s.series_id)
       if (s.dim4 == 1020) and ('Run4_cmrr_mbep2d_bold_25mm_GRAPPA3_PF78_ME3_MB2' in s.series_description):
            info[func_task07].append(s.series_id)
       if (s.dim4 == 15) and ('Run4_inv_cmrr_mbep2d_bold_25mm_GRAPPA3_PF78_ME3_MB2' in s.series_description):
            info[func_task08].append(s.series_id)
       if (s.dim4 > 171) and ('Run5_ep2d_bold_ptx_20220214_VERSE_5pulses_8ms_FA40' in s.series_description):
            info[func_task09].append(s.series_id)
       if (s.dim4 == 5) and ('Run5_inv_ep2d_bold_ptx_20220214_VERSE_5pulses_8ms_FA40' in s.series_description):
            info[func_task10].append(s.series_id)

    return info
