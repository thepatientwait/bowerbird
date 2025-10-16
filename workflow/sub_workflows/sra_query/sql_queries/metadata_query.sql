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
    geo_loc_name_country_calc,
    geo_loc_name_country_continent_calc,
    geo_loc_name_sam,
    organism,
    libraryselection,
    librarysource,
    sra_study,
    avgspotlen
FROM
    `nih-sra-datastore.sra.metadata`
WHERE
    (librarysource = 'METAGENOMIC'
        OR organism IN (SELECT sci_name FROM metadata.emp))
    AND librarysource NOT IN ('VIRAL RNA', 'METATRANSCRIPTOMIC', 'TRANSCRIPTOMIC')
    AND consent = 'public'
    {start_date_filter}
    {target_filter};