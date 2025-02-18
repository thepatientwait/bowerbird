-- Loop counter
DECLARE counter int64 default 1;

-- Create intermediate table with metagenome taxonomy entries
CREATE OR REPLACE TABLE 
    metadata.emp AS
WITH 
    cte AS (
        SELECT 
            1 AS xlevel, tax_id, parent_id, sci_name
        FROM 
            `nih-sra-datastore.sra_tax_analysis_tool.taxonomy` e
        WHERE 
            e.sci_name = 'metagenomes'
        )
SELECT *
FROM cte;

WHILE EXISTS (
    SELECT c.*
    FROM 
        metadata.emp p
    INNER JOIN 
        `nih-sra-datastore.sra_tax_analysis_tool.taxonomy` c 
        ON p.tax_id = c.parent_id
    WHERE 
        p.xlevel = counter
    )
DO
-- Insert next level
INSERT INTO 
    metadata.emp ( xlevel, tax_id, parent_id, sci_name )
SELECT 
    counter + 1 AS xlevel, c.tax_id, c.parent_id, c.sci_name
FROM 
    metadata.emp p
INNER JOIN 
    `nih-sra-datastore.sra_tax_analysis_tool.taxonomy` c ON p.tax_id = c.parent_id
WHERE 
    p.xlevel = counter;

SET counter = counter + 1;

-- Loop safely
IF counter > 10
THEN
    BREAK;
END IF;
END WHILE;

-- Select for download
SELECT *
FROM metadata.emp;