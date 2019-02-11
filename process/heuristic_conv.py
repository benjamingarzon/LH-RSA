
import os

def create_key(template, outtype=('nii.gz','dicom'), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return (template, outtype, annotation_classes)


def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    allowed template fields - follow python string module:
    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    t1 = create_key('anat/sub-{subject}_T1w')
    sequence = create_key('func/sub-{subject}_task-sequence_run-{item:02d}_bold')    
    fmap = create_key('fmap/sub-{subject}_epi')
    
    #info = {t1:[], sequence:[], fmap:[]}
    info = {t1:[]}
    
    for idx, s in enumerate(seqinfo):
        print(idx, s)
        if (s.dim3 == 225) and (s.dim4 == 1) and ('T1 Volume' in s.protocol_name):
            info[t1] = [s.series_id]
#        if (s.dim3 == 22440) and ('1101' in s.series_id):
#            info[sequence].append({'item': s.series_id})
#        if (s.dim3 == 22440) and ('1201' in s.series_id):
#            info[sequence].append({'item': s.series_id})
#        if (s.dim3 == 22440) and ('1301' in s.series_id):
#            info[sequence].append({'item': s.series_id})
#        if (s.dim3 == 22440) and ('1401' in s.series_id):
#            info[sequence].append({'item': s.series_id})
#        if ('501' in s.series_id):
#            info[fmap] = [s.series_id]

    return info
