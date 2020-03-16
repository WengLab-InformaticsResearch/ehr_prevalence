-- Export condition, drug, procedure, and demographics data from SQL Server Management Studio for EHR prevalence analysis
--
-- Setup:
-- Set file output format to tab delimited:
--     In SQL Server Management Studio: Tools > Options > Query Results > SQL Server > Results to Text >
-- 	       Output format: tab delimited
--         Include column headers in the result set: enabled
--     Restart SSMS for new settings to take effect


-- Prevent the count from showing up in the text file results
SET NOCOUNT ON;


-- Load the iatrogenic codes into temporary tables: #iatrogenic_codes and #iatrogenic_codes_with_desc
-- Iatrogenic codes and their descendants will be excluded from the analysis
:r D:\cohd\ehr_prevalence\load_iatrogenic_tsql.sql


-- -------------------------------------------------------------
-- Export conditions, drugs, and procedures, excluding iatrogenic codes
-- -------------------------------------------------------------

-- Export person ID, start date, and concept IDs for conditions, drugs, and procedures
:OUT D:\cohd\data\unique_patient_concept_pairs_date.txt
SELECT DISTINCT co.person_id, YEAR(co.condition_start_date) AS date, co.condition_concept_id AS concept_id
FROM condition_occurrence co
JOIN concept c ON co.condition_concept_id = c.concept_id
LEFT JOIN #iatrogenic_codes_with_desc i ON c.concept_id = i.concept_id
WHERE condition_concept_id != 0 
	AND c.domain_id = 'Condition'	-- Make sure we only get conditions from the condition_occurrence table
	AND i.concept_id IS NULL		-- Make sure condition is not an iatrogenic code
UNION ALL
SELECT DISTINCT de.person_id, YEAR(de.drug_exposure_start_date) AS date, de.drug_concept_id AS concept_id
FROM drug_exposure de
JOIN concept c ON de.drug_concept_id = c.concept_id
LEFT JOIN #iatrogenic_codes_with_desc i ON c.concept_id = i.concept_id
WHERE drug_concept_id != 0 
	AND c.domain_id = 'Drug'	-- Make sure we only get conditions from the condition_occurrence table
	AND i.concept_id IS NULL	-- Make sure condition is not an iatrogenic code
UNION ALL
SELECT DISTINCT po.person_id, YEAR(po.procedure_date) AS date, po.procedure_concept_id AS concept_id
FROM procedure_occurrence po
JOIN concept c ON po.procedure_concept_id = c.concept_id
LEFT JOIN #iatrogenic_codes_with_desc i ON c.concept_id = i.concept_id
WHERE procedure_concept_id != 0 
	AND c.domain_id = 'Procedure'	-- Make sure we only get conditions from the condition_occurrence table
	AND i.concept_id IS NULL		-- Make sure condition is not an iatrogenic code;

-- Export demographics from the person table
:OUT D:\cohd\data\person.txt
SELECT person_id, gender_concept_id, race_concept_id, ethnicity_concept_id
FROM person;

-- -------------------------------------------------------------
-- Concept hierarchy
-- Get observed concepts and their ancestors from the desired vocabularies
-- -------------------------------------------------------------
-- Get observed concepts from conditions, drugs, and procedures
SELECT DISTINCT condition_concept_id AS concept_id
INTO #observed_condition_concepts
FROM condition_occurrence o
JOIN concept c ON o.condition_concept_id = c.concept_id
LEFT JOIN #iatrogenic_codes_with_desc i ON c.concept_id = i.concept_id
WHERE c.domain_id = 'Condition' AND i.concept_id IS NULL;

SELECT DISTINCT drug_concept_id AS concept_id
INTO #observed_drug_concepts
FROM drug_exposure o
JOIN concept c ON o.drug_concept_id = c.concept_id
LEFT JOIN #iatrogenic_codes_with_desc i ON c.concept_id = i.concept_id
WHERE c.domain_id = 'Drug' AND i.concept_id IS NULL;

SELECT DISTINCT procedure_concept_id AS concept_id
INTO #observed_procedure_concepts
FROM procedure_occurrence o
JOIN concept c ON o.procedure_concept_id = c.concept_id
LEFT JOIN #iatrogenic_codes_with_desc i ON c.concept_id = i.concept_id
WHERE c.domain_id = 'Procedure' AND i.concept_id IS NULL;

-- Get condition concepts and their ancestors
SELECT * 
INTO #hierarchical_condition_concepts
FROM
	((SELECT *
	FROM #observed_condition_concepts)
	UNION
	(SELECT DISTINCT ca.ancestor_concept_id
	FROM #observed_condition_concepts x
	JOIN concept_ancestor ca ON ca.descendant_concept_id = x.concept_id
	JOIN concept c ON ca.ancestor_concept_id = c.concept_id
	WHERE c.domain_id = 'Condition' AND c.vocabulary_id = 'SNOMED')
	UNION
	(SELECT DISTINCT ca.ancestor_concept_id
	FROM #observed_condition_concepts x
	JOIN concept_ancestor ca ON ca.descendant_concept_id = x.concept_id
	JOIN concept c ON ca.ancestor_concept_id = c.concept_id
	LEFT JOIN concept_relationship cr ON (cr.concept_id_1 = c.concept_id AND cr.relationship_id = 'MedDRA - SNOMED eq')
	WHERE c.domain_id = 'Condition' AND c.vocabulary_id = 'MedDRA' AND cr.relationship_id IS NULL)) y
;


-- Get drug concepts and their ancestors
SELECT * 
INTO #hierarchical_drug_concepts
FROM
	((SELECT *
	FROM #observed_drug_concepts)
	UNION
	(SELECT DISTINCT ca.ancestor_concept_id
	FROM #observed_drug_concepts x
	JOIN concept_ancestor ca ON ca.descendant_concept_id = x.concept_id
	JOIN concept c ON ca.ancestor_concept_id = c.concept_id
	WHERE c.domain_id = 'Drug' AND c.vocabulary_id = 'RxNorm' AND c.concept_class_id IN ('Ingredient', 'Clinical Drug Form', 'Clinical Drug Comp', 'Clinical Drug'))
	UNION
	(SELECT DISTINCT ca.ancestor_concept_id
	FROM #observed_drug_concepts x
	JOIN concept_ancestor ca ON ca.descendant_concept_id = x.concept_id
	JOIN concept c ON ca.ancestor_concept_id = c.concept_id
	LEFT JOIN concept_relationship cr ON (cr.concept_id_1 = c.concept_id AND cr.relationship_id IN ('ATC - RxNorm', 'ATC - RxNorm name'))
	WHERE c.domain_id = 'Drug' AND c.vocabulary_id = 'ATC' AND cr.relationship_id IS NULL)) y
;

-- Get procedure concepts and their ancestors
SELECT * 
INTO #hierarchical_procedure_concepts
FROM
	((SELECT *
	FROM #observed_procedure_concepts)
	UNION
	(SELECT DISTINCT ca.ancestor_concept_id
	FROM #observed_procedure_concepts x
	JOIN concept_ancestor ca ON ca.descendant_concept_id = x.concept_id
	JOIN concept c ON ca.ancestor_concept_id = c.concept_id
	WHERE c.domain_id = 'Procedure' AND c.vocabulary_id IN ('SNOMED', 'ICD10PCS'))
	UNION
	(SELECT DISTINCT ca.ancestor_concept_id
	FROM #observed_procedure_concepts x
	JOIN concept_ancestor ca ON ca.descendant_concept_id = x.concept_id
	JOIN concept c ON ca.ancestor_concept_id = c.concept_id
	LEFT JOIN concept_relationship cr ON (cr.concept_id_1 = c.concept_id AND cr.relationship_id = 'MedDRA - SNOMED eq')
	WHERE c.domain_id = 'Procedure' AND c.vocabulary_id = 'MedDRA' AND cr.relationship_id IS NULL)) y
;

-- Get the vocabularies used in the condition/drug/procedure tables and for their hierarchy
declare @condition_vocabs TABLE (vocabulary_id varchar(20));
declare @drug_vocabs TABLE (vocabulary_id varchar(20));
declare @procedure_vocabs TABLE (vocabulary_id varchar(20));
INSERT INTO @condition_vocabs (vocabulary_id)
	(SELECT * 
	FROM 
	-- Vocabularies used in hierarchy
	(values ('SNOMED'), ('MedDRA')) v(vocabulary_id))
	UNION
	-- Vocabularies in observational tables
	(SELECT DISTINCT vocabulary_id
	FROM condition_occurrence co
	JOIN concept c ON co.condition_concept_id = c.concept_id
	WHERE c.domain_id = 'Condition');
INSERT INTO @drug_vocabs (vocabulary_id)
	(SELECT * 
	FROM 
	-- Vocabularies used in hierarchy
	(values ('RxNorm'), ('ATC')) v(vocabulary_id))
	UNION
	-- Vocabularies in observational tables
	(SELECT DISTINCT vocabulary_id
	FROM drug_exposure de2
	JOIN concept c2 ON de2.drug_concept_id = c2.concept_id
	WHERE c2.domain_id = 'Drug');
INSERT INTO @procedure_vocabs (vocabulary_id)
	(SELECT * 
	FROM 
	-- Vocabularies used in hierarchy
	(values ('SNOMED'), ('MedDRA'), ('ICD10PCS')) v(vocabulary_id))
	UNION
	-- Vocabularies in observational tables
	(SELECT DISTINCT vocabulary_id
	FROM procedure_occurrence po
	JOIN concept c ON po.procedure_concept_id = c.concept_id
	WHERE c.domain_id = 'Procedure');

-- Export the hierarchical concepts and their direct descendants from the desired vocabularies
-- Note: some concepts can be descendants with 0 levels of separation from the ancestor
-- Note 2: Using variables to store the vocabs instead of using a subquery because the 
-- subquery was extremely slow despite being non-correlated. 
:OUT D:\cohd\data\concept_descendants_direct.txt
SELECT DISTINCT concept_id, descendant_concept_id
FROM
	((SELECT ca.ancestor_concept_id AS concept_id, ca.descendant_concept_id
	FROM #hierarchical_condition_concepts hc
	JOIN concept_ancestor ca ON ca.ancestor_concept_id = hc.concept_id
	JOIN concept c ON ca.descendant_concept_id = c.concept_id
	-- Make sure the descendant concept is one of the observed concepts or its ancestors
	JOIN #hierarchical_condition_concepts hc2 ON ca.descendant_concept_id = hc2.concept_id 
	WHERE 
		-- Don't include the self-ancestor relationship
		ca.ancestor_concept_id != ca.descendant_concept_id
		-- Want only immediate descendants (or 0 level descendants)
		AND ca.min_levels_of_separation <= 1 
		-- Want only descendants from the same domain
		AND c.domain_id = 'Condition'
		-- Want only descendants from the specified vocabs
		AND c.vocabulary_id IN 
			(SELECT * FROM @condition_vocabs))
	UNION
	(SELECT ca.ancestor_concept_id AS concept_id, ca.descendant_concept_id
	FROM #hierarchical_drug_concepts hc
	JOIN concept_ancestor ca ON ca.ancestor_concept_id = hc.concept_id
	JOIN concept c ON ca.descendant_concept_id = c.concept_id
	-- Make sure the descendant concept is one of the observed concepts or its ancestors
	JOIN #hierarchical_drug_concepts hc2 ON ca.descendant_concept_id = hc2.concept_id
	WHERE 
		-- Don't include the self-ancestor relationship
		ca.ancestor_concept_id != ca.descendant_concept_id
		-- Want only immediate descendants (or 0 level descendants)
		AND ca.min_levels_of_separation <= 1 
		-- Want only descendants from the same domain
		AND c.domain_id = 'Drug'
		-- Want only descendants from the specified vocabs
		AND c.vocabulary_id IN
			(SELECT * FROM @drug_vocabs))
	UNION
	(SELECT ca.ancestor_concept_id AS concept_id, ca.descendant_concept_id
	FROM #hierarchical_procedure_concepts hc
	JOIN concept_ancestor ca ON ca.ancestor_concept_id = hc.concept_id
	JOIN concept c ON ca.descendant_concept_id = c.concept_id
	-- Make sure the descendant concept is one of the observed concepts or its ancestors
	JOIN #hierarchical_procedure_concepts hc2 ON ca.descendant_concept_id = hc2.concept_id
	WHERE 
		-- Don't include the self-ancestor relationship
		ca.ancestor_concept_id != ca.descendant_concept_id
		-- Want only immediate descendants (or 0 level descendants)
		AND ca.min_levels_of_separation <= 1 
		-- Want only descendants from the same domain
		AND c.domain_id = 'Procedure'
		-- Want only descendants from the specified vocabs
		AND c.vocabulary_id IN
			(SELECT * FROM @procedure_vocabs))) tmp
ORDER BY concept_id, descendant_concept_id
;

-- Export all observed descendants of each each hierarchical concept
-- Note: some concepts can be descendants with 0 levels of separation from the ancestor
-- Note 2: Using variables to store the vocabs instead of using a subquery because the
-- subquery was extremely slow despite being non-correlated.
:OUT D:\cohd\data\concept_descendants_all_observed.txt
SELECT concept_id, descendant_concept_id
FROM
        ((SELECT ca.ancestor_concept_id AS concept_id, ca.descendant_concept_id
        FROM #hierarchical_condition_concepts hc
        JOIN concept_ancestor ca ON ca.ancestor_concept_id = hc.concept_id
        JOIN concept c ON ca.descendant_concept_id = c.concept_id
        -- Make sure the descendant concept is one of the observed concepts
        JOIN #observed_condition_concepts oc ON ca.descendant_concept_id = oc.concept_id
        WHERE
                -- Don't include the self-ancestor relationship
                ca.ancestor_concept_id != ca.descendant_concept_id
                -- Want only descendants from the same domain
                AND c.domain_id = 'Condition'
                -- Want only descendants from the specified vocabs
                AND c.vocabulary_id IN
                        (SELECT * FROM @condition_vocabs))
        UNION
        (SELECT ca.ancestor_concept_id AS concept_id, ca.descendant_concept_id
        FROM #hierarchical_drug_concepts hc
        JOIN concept_ancestor ca ON ca.ancestor_concept_id = hc.concept_id
        JOIN concept c ON ca.descendant_concept_id = c.concept_id
        -- Make sure the descendant concept is one of the observed concepts
        JOIN #observed_drug_concepts oc ON ca.descendant_concept_id = oc.concept_id
        WHERE
                -- Don't include the self-ancestor relationship
                ca.ancestor_concept_id != ca.descendant_concept_id
                -- Want only descendants from the same domain
                AND c.domain_id = 'Drug'
                -- Want only descendants from the specified vocabs
                AND c.vocabulary_id IN
                        (SELECT * FROM @drug_vocabs))
        UNION
        (SELECT ca.ancestor_concept_id AS concept_id, ca.descendant_concept_id
        FROM #hierarchical_procedure_concepts hc
        JOIN concept_ancestor ca ON ca.ancestor_concept_id = hc.concept_id
        JOIN concept c ON ca.descendant_concept_id = c.concept_id
        -- Make sure the descendant concept is one of the observed concepts
        JOIN #observed_procedure_concepts oc ON ca.descendant_concept_id = oc.concept_id
        WHERE
                -- Don't include the self-ancestor relationship
                ca.ancestor_concept_id != ca.descendant_concept_id
                -- Want only descendants from the same domain
                AND c.domain_id = 'Procedure'
                -- Want only descendants from the specified vocabs
                AND c.vocabulary_id IN
                        (SELECT * FROM @procedure_vocabs))) tmp
ORDER BY concept_id, descendant_concept_id
;


-- Export the concept definitions
-- Note: use union instead of union all because 0 is in each domain
:OUT D:\cohd\data\concepts.txt
SELECT *
FROM
	(SELECT concept_id FROM #hierarchical_condition_concepts
	UNION
	SELECT concept_id FROM #hierarchical_drug_concepts
	UNION
	SELECT concept_id FROM #hierarchical_procedure_concepts
	UNION
	SELECT DISTINCT gender_concept_id AS concept_id FROM person
	UNION
	SELECT DISTINCT race_concept_id AS concept_id FROM person
	UNION
	SELECT DISTINCT ethnicity_concept_id AS concept_id FROM person) concept_ids
LEFT JOIN concept ON concept_ids.concept_id = concept.concept_id
ORDER BY concept.concept_id ASC;


-- Return to normal settings
SET NOCOUNT OFF;

