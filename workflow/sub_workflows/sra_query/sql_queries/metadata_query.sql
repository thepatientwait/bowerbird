SELECT 
    acc,
    experiment,
    biosample,
    bioproject,
    mbases,
    mbytes,
    assay_type,
    releasedate,
    collection_date_sam,
    organism,
    libraryselection,
    librarysource,
    sra_study,
    avgspotlen,
    attributes
FROM
    `nih-sra-datastore.sra.metadata`
WHERE
    (librarysource = 'METAGENOMIC'
        OR organism IN (SELECT sci_name FROM metadata.emp))
    AND librarysource NOT IN ('VIRAL RNA', 'METATRANSCRIPTOMIC', 'TRANSCRIPTOMIC')
    AND consent = 'public'
    {start_date_filter}
    {target_filter};