#!/bin/bash
# ---- Usage Check ----
if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <sample_file> <pos_file> <H5_FILE_NAME>"
  exit 1
fi

SAMPLE_FILE="$1"

POS_FILE="$2"

H5_FILE_NAME="$3"

if [[ ! -f "$SAMPLE_FILE" ]]; then
  echo "❌ Error: File '$SAMPLE_FILE' not found."
  exit 1
fi

if [[ ! -f "$POS_FILE" ]]; then
  echo "❌ Error: File '$POS_FILE' not found."
  exit 1
fi

# ---- Configuration ----
DB_NAME="<INSERT DB_NAME>"
DB_USER="<INSERT USERNAME>"
DB_HOST="<INSERT URL>"
DB_PORT="<INSERT PORT>"
DB_PASSWORD="<INSERT PASSWORD>"
CSV_DIR="/path/to/csv/files"
LOG_FILE="upload_log.txt"

# ---- Start Upload ----
echo "Starting bulk upload at $(date)" > "$LOG_FILE"

CVTERM_WHOLE_PLANT='whole plant'
CV__PLANT_ANATOMY='plant_anatomy'
ORGANISM_NAME='Japonica nipponbare'
ORGANISM_GENUS='Oryza'
ORGANISM_SPECIES='Oryza sativa'
DB='1k1'



VARIANTSET_NAME='1k1'
CVTERM_SNP='SNP'
CVTERM_CHR='chromosome'
VARIANT_TYPE_CV='sequence'  # or use 'SO' or your CV name

OUTPUT_STOCK_SAMPLE_CSV="StockSample.csv"
OUTPUT_SAMPLE_VARIETYSET="SampleVarietySet.csv"
OUTPUT_SNP_FEATURE="SnpFeature.csv"
OUTPUT_VARIANT_VARIANTSET="Variant_variantset.csv"
OUTPUT_SNP_FEATURELOC="Snp_FeatureLoc.csv"


# ---- CVTERM Lookup or Insert (CV = Plant anatomy and CVterm ----
echo "Looking up cvterm_id for '$CVTERM_WHOLE_PLANT' from '$CV__PLANT_ANATOMY'..." | tee -a "$LOG_FILE"

CVTERM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
SELECT cvterm_id FROM cvterm JOIN cv ON cvterm.cv_id = cv.cv_id 
WHERE cvterm.name = '$CVTERM_WHOLE_PLANT' AND cv.name = '$CV__PLANT_ANATOMY' LIMIT 1;
")

CVTERM_ID=$(echo $CVTERM_ID | xargs)

if [[ -z "$CVTERM_ID" ]]; then
  echo "⚠️ cvterm_id not found, inserting cv and cvterm..." | tee -a "$LOG_FILE"

  # Ensure CV exists
  CV_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
  INSERT INTO cv(name) 
  VALUES ('$CV__PLANT_ANATOMY') 
  ON CONFLICT (name) DO UPDATE SET name=EXCLUDED.name 
  RETURNING cv_id;
  ")

  CV_ID=$(echo $CV_ID | xargs)

  # Insert cvterm
  CVTERM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
  INSERT INTO cvterm(name, cv_id, is_obsolete) 
  VALUES ('$CVTERM_WHOLE_PLANT', $CV_ID, false)
  RETURNING cvterm_id;
  ")
fi

CVTERM_ID=$(echo $CVTERM_ID | xargs)
echo "✅ GOT CVTERM_ID: $CVTERM_ID" | tee -a "$LOG_FILE"

# ---- Organism Lookup or Insert ----
echo "Looking up organism_id for '$ORGANISM_NAME'..." | tee -a "$LOG_FILE"

ORGANISM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
SELECT organism_id FROM organism WHERE common_name = '$ORGANISM_NAME' LIMIT 1;
")

if [[ -z "$ORGANISM_ID" ]]; then
  echo "⚠️ organism_id not found, inserting organism..." | tee -a "$LOG_FILE"

  ORGANISM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
  INSERT INTO organism (genus, species, common_name)
  VALUES ('$ORGANISM_GENUS', '$ORGANISM_SPECIES', '$ORGANISM_NAME')
  RETURNING organism_id;
  ")
fi

ORGANISM_ID=$(echo $ORGANISM_ID | xargs)
echo "✅ GOT ORGANISM_ID: $ORGANISM_ID" | tee -a "$LOG_FILE"

# FEATURE TERM
echo "Looking up FEATURE CVTERM  for CV = '$VARIANT_TYPE_CV'...and cvterm: '$CVTERM_CHR'" | tee -a "$LOG_FILE"
FEATURE_CVTERM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
	SELECT cvterm_id 
	FROM cvterm 
	JOIN cv ON cv.cv_id = cvterm.cv_id 
	WHERE cv.name = '$VARIANT_TYPE_CV' AND cvterm.name = '$CVTERM_CHR'
	LIMIT 1;
	")

FEATURE_CVTERM_ID=$(echo $FEATURE_CVTERM_ID | xargs)

# ---- DB Lookup or Insert ----
echo "Looking up db_id for '$DB'..." | tee -a "$LOG_FILE"

DB_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
SELECT db_id FROM db WHERE name = '$DB' LIMIT 1;
")

if [[ -z "$DB_ID" ]]; then
  echo "⚠️ db_id not found, inserting db..." | tee -a "$LOG_FILE"

  DB_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
  INSERT INTO db(name)
  VALUES ('$DB')
  RETURNING db_id;
  ")
fi

DB_ID=$(echo $DB_ID | xargs)
echo "✅ GOT DB_ID: $DB_ID" | tee -a "$LOG_FILE"

# ---- VariantSet Lookup or Insert ----

echo "Looking up variantset_id for variantset name '$VARIANTSET_NAME'..." | tee -a "$LOG_FILE"

VARIANTSET_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
SELECT variantset_id FROM variantset WHERE name = '$VARIANTSET_NAME' LIMIT 1;")

if [[ -z "$VARIANTSET_ID" ]]; then
  echo "⚠️ variantset_id not found." | tee -a "$LOG_FILE"
  echo "Looking up cvterm_id for variant type '$CVTERM_SNP' in CV '$VARIANT_TYPE_CV'..." | tee -a "$LOG_FILE"

	VARIANTTYPE_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
	SELECT cvterm_id 
	FROM cvterm 
	JOIN cv ON cv.cv_id = cvterm.cv_id 
	WHERE cv.name = '$VARIANT_TYPE_CV' AND cvterm.name = '$CVTERM_SNP'
	LIMIT 1;
	")

	VARIANTTYPE_ID=$(echo $VARIANTTYPE_ID | xargs)

	if [[ -z "$VARIANTTYPE_ID" ]]; then
	  echo "❌ Could not find cvterm_id for variant type '$CVTERM_SNP' in '$VARIANT_TYPE_CV'" | tee -a "$LOG_FILE"
	  exit 1
	fi

  echo "✅ GOT VARIANTTYPE_ID: $VARIANTTYPE_ID" | tee -a "$LOG_FILE"

  echo "Adding variantset_id." | tee -a "$LOG_FILE"
  
  echo "Fetching current max variant_set_id..." | tee -a "$LOG_FILE"

  MAX_VARIANT_SET_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -t -A -c \
	  "SELECT COALESCE(MAX(variantset_id), 0) FROM variantset;")

	MAX_VARIANT_SET_ID=$(echo "$MAX_VARIANT_SET_ID" | xargs)
	NEXT_VARIANT_SET_ID=$((MAX_VARIANT_SET_ID + 1))

	PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -c \
	  "SELECT setval('variantset_variantset_id_seq', $NEXT_VARIANT_SET_ID, false);"


  VARIANTSET_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
  INSERT INTO variantset(name, variant_type_id)
  VALUES ('$VARIANTSET_NAME', $VARIANTTYPE_ID)
  RETURNING variantset_id;
  ")
fi

VARIANTSET_ID=$(echo $VARIANTSET_ID | xargs)
echo "✅ GOT VARIANTSET_ID: $VARIANTSET_ID" | tee -a "$LOG_FILE"


# ---- Add Sample to database ----

echo "Starting sample loading..." | tee "$LOG_FILE"



MAX_STOCK_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -c "SELECT COALESCE(MAX(stock_id), 0) FROM stock;" | xargs)

echo "🔢 Max stock: $MAX_STOCK_ID"

	# Update the sequence (replace sequence name with actual name)
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('stock_stock_id_seq', $MAX_STOCK_ID + 1, false);"

STOCK_ID=$((MAX_STOCK_ID +1))

echo "✅ (STOCK_SAMPLE) Sequence updated to start at $((STOCK_ID))"



MAX_STOCK_SAMPLE_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -c "SELECT COALESCE(MAX(stock_sample_id), 0) FROM stock_sample;" | xargs)

echo "🔢 Max stock_sample_id: $MAX_STOCK_SAMPLE_ID"

	# Update the sequence (replace sequence name with actual name)
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('stock_sample_stock_sample_id_seq', $MAX_STOCK_SAMPLE_ID + 1, false);"

STOCK_SAMPLE_ID=$((MAX_STOCK_SAMPLE_ID +1))

echo "✅ (STOCK_SAMPLE) Sequence updated to start at $((STOCK_SAMPLE_ID))"


echo "stock_sample_id, stock_id,dbxref_id,hdf5_index,tmp_oldstock_id " > "$OUTPUT_STOCK_SAMPLE_CSV"

# --- Get max sample_varietyset_id ---
MAX_SAMPLE_VARIETYSET_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -c "SELECT COALESCE(MAX(sample_varietyset_id), 0) FROM sample_varietyset;" | xargs)

echo "🔢 Max sample_varietyset_id: $MAX_SAMPLE_VARIETYSET_ID"

# --- Update the sequence ---
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('sample_run_sample_run_id_seq', $((MAX_SAMPLE_VARIETYSET_ID + 1)), false);
"

SAMPLE_VARIETYSET_ID=$((MAX_SAMPLE_VARIETYSET_ID +1))

echo "✅ (SAMPLE VARIETY_SET)  Sequence updated to start at $((SAMPLE_VARIETYSET_ID))"


echo "sample_varietyset_id, stock_sample_id, db_id,hdf5_index" > "$OUTPUT_SAMPLE_VARIETYSET"


COUNTER=0
HDF_COUNTER=0

while IFS= read -r LINE || [ -n "$LINE" ]; do
  echo "Processing line: $LINE" | tee -a "$LOG_FILE"

     # ---- Add / Check for Stock Record ----

  STOCK_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
    SELECT stock_id FROM stock
    WHERE uniquename = '$LINE'
	  AND name = '$LINE'
      AND type_id = $CVTERM_ID
      AND organism_id = $ORGANISM_ID
    LIMIT 1;
  " | xargs)

  if [[ -z "$STOCK_ID" ]]; then
    echo "⚠️ Stock '$LINE' not found. Checking Stock Prop" | tee -a "$LOG_FILE"

    # check in Stock_property table
    STOCK_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
    select sp.stock_id  from stockprop sp where sp.value like '$LINE' 
    LIMIT 1;
    " | xargs)
    
    if [[ -z "$STOCK_ID" ]]; then
      STOCK_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
        INSERT INTO stock (name, uniquename, type_id, organism_id)
        VALUES ('$LINE', '$LINE', $CVTERM_ID, $ORGANISM_ID)
        RETURNING stock_id; " | xargs)

        echo "✅ Inserted new stock_id: $STOCK_ID at $HDF_COUNTER" | tee -a "$LOG_FILE"
    else
      echo "✅ Found existing stock_id in stockprop: $STOCK_ID at $HDF_COUNTER" | tee -a "$LOG_FILE"
    fi
  fi

  # --- DBXREF INSERT OR SELECT ---
	echo "Looking up dbxref_id for db_id=$DB_ID and accession='$LINE'..." | tee -a "$LOG_FILE"

	DBXREF_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
	  SELECT dbxref_id FROM dbxref
	  WHERE db_id = $DB_ID AND accession = '$LINE'
	  LIMIT 1;
	" | xargs)

	if [[ -z "$DBXREF_ID" ]]; then
	  echo "⚠️ dbxref not found. Inserting..." | tee -a "$LOG_FILE"

	  DBXREF_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
		INSERT INTO dbxref (db_id, accession, version)
		VALUES ($DB_ID, '$LINE', 1)
		RETURNING dbxref_id;
	  " | xargs)

	  echo "✅ Inserted new dbxref_id: $DBXREF_ID" | tee -a "$LOG_FILE"
	else
	  echo "✅ Found existing dbxref_id: $DBXREF_ID" | tee -a "$LOG_FILE"
	fi
	
  
  echo "$STOCK_SAMPLE_ID, $STOCK_ID,$DBXREF_ID, $HDF_COUNTER, $HDF_COUNTER " >> "$OUTPUT_STOCK_SAMPLE_CSV"
  
  echo "$SAMPLE_VARIETYSET_ID, $STOCK_SAMPLE_ID,$DB_ID, $HDF_COUNTER" >> "$OUTPUT_SAMPLE_VARIETYSET"
  
  
  ((STOCK_SAMPLE_ID++))
   ((SAMPLE_VARIETYSET_ID++))
  ((HDF_COUNTER++))
 

done < "$SAMPLE_FILE"


# ---- SNP Feature and VariantSet Upload ----

COUNTER=0

MAX_SNP_FEATURE_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -c "SELECT COALESCE(MAX(snp_feature_id), 0) FROM snp_feature;" | xargs)

echo "🔢 Max stock: $MAX_SNP_FEATURE_ID"

	# Update the sequence (replace sequence name with actual name)
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('snp_feature_snp_feature_id_seq', $MAX_SNP_FEATURE_ID + 1, false);"

SNP_FEATURE_ID=$((MAX_STMAX_SNP_FEATURE_IDOCK_ID +1))

echo "✅ (SNP_FEATURE) Sequence updated to start at $((SNP_FEATURE_ID))"

MAX_SNP_FEATURE_LOC_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -c "SELECT COALESCE(MAX(snp_featureloc_id), 0) FROM snp_featureloc;" | xargs)

echo "🔢 SNP_FEATURELOC : $MAX_SNP_FEATURE_LOC_ID"

	# Update the sequence (replace sequence name with actual name)
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('snp_featureloc_snp_featureloc_id_seq', $MAX_SNP_FEATURE_LOC_ID + 1, false);"

SNP_FEATURELOC_ID=$((MAX_SNP_FEATURE_LOC_ID +1))

echo "✅ (SNP_FEATURE_LOC) Sequence updated to start at $((SNP_FEATURELOC_ID))"


MAX_VARIANT_VARIANTSET_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -c "SELECT COALESCE(MAX(variant_variantset_id), 0) FROM variant_variantset;" | xargs)

echo "🔢 VARIANT_VARIANTSET_ID : $MAX_VARIANT_VARIANTSET_ID"

	# Update the sequence (replace sequence name with actual name)
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('variant_variantset_variant_variantset_id_seq', $MAX_VARIANT_VARIANTSET_ID + 1, false);"

VARIANT_VARIANTSET_ID=$((MAX_VARIANT_VARIANTSET_ID +1))

echo "✅ (VARIANT_VARIANTSET) Sequence updated to start at $((VARIANT_VARIANTSET_ID))"




echo "snp_feature_id, variantset_id" > $OUTPUT_SNP_FEATURE

echo "variant_variantset_id, variant_feature_id, variantset_id, hdf5_index" > $OUTPUT_VARIANT_VARIANTSET

echo "snp_featureloc_id, organism_id, srcfeature_id, snp_feature_id, position, refcall" > $OUTPUT_SNP_FEATURELOC

while IFS=$'\t' read -r chrom posit refc rest; do
	if [[ "$chrom" != "$POS" ]]; then
		NAME="Chr$chrom"
		UNIQUENAME="chr0$chrom"
		
		echo "Inserting new feature for $NAME $UNIQUENAME $COUNTER..."
		
		echo "NAME: $NAME"
		echo "UNIQUENAME: $UNIQUENAME"
		echo "FEATURE_CVTERM_ID: $FEATURE_CVTERM_ID"
		echo "ORGANISM_ID: $ORGANISM_ID"

		
		 FEATURE_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
			  SELECT feature_id 
			  FROM feature
			  WHERE name = '$NAME' 
				AND uniquename = '$UNIQUENAME'
				AND type_id = $FEATURE_CVTERM_ID
				AND organism_id = $ORGANISM_ID
			  LIMIT 1;
			")

			FEATURE_ID=$(echo "$FEATURE_ID" | xargs)
	
		POS=$chrom
	fi
	
	echo "$SNP_FEATURE_ID, $VARIANTSET_ID" >> $OUTPUT_SNP_FEATURE
	
	echo "$VARIANT_VARIANTSET_ID, $SNP_FEATURE_ID,  $VARIANTSET_ID, $COUNTER" >> $OUTPUT_VARIANT_VARIANTSET
	
	FINAL_POS=$((posit - 1))
	
	echo "$SNP_FEATURELOC_ID, $ORGANISM_ID,  $FEATURE_ID, $SNP_FEATURE_ID, $FINAL_POS, $refc " >> $OUTPUT_SNP_FEATURELOC
	
	((SNP_FEATURE_ID++))
	((VARIANT_VARIANTSET_ID++))
	((SNP_FEATURELOC_ID++))
	#echo "adding $counter"
	 

	((COUNTER++))
done < "$POS_FILE"

PLATFORM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
    select platform_id from platform p where p.variantset_id = $VARIANTSET_ID and p.db_id =$DB_ID
    LIMIT 1;
  " | xargs)

  if [[ -z "$PLATFORM_ID" ]]; then
    echo "⚠️ PLATFORM_ID not found. Inserting... Db=$DB_ID and VariantSetId=$VARIANTSET_ID" | tee -a "$LOG_FILE"

    PLATFORM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
      INSERT INTO platform (variantset_id, db_id)
      VALUES ($VARIANTSET_ID, $DB_ID)
      RETURNING platform_id;
    " | xargs)

    echo "✅ Inserted new PLATFORM_ID: $PLATFORM_ID" | tee -a "$LOG_FILE"
  else
    echo "✅ Found existing PLATFORM_ID: $PLATFORM_ID" | tee -a "$LOG_FILE"
  fi


  GENOTYPE_RUN_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
    select gr.genotype_run_id from genotype_run gr where gr.platform_id =PLATFORM_ID 
    LIMIT 1;
  " | xargs)

  if [[ -z "$GENOTYPE_RUN_ID" ]]; then
    echo "⚠️ PLATFORM_ID not found. Inserting... Db=$DB_ID and VariantSetId=$VARIANTSET_ID" | tee -a "$LOG_FILE"

    GENOTYPE_RUN_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -c "
      INSERT INTO platform (platform_id)
      VALUES ($PLATFORM_ID, '$H5_FILE_NAME')
      RETURNING genotype_run_id;
    " | xargs)

    echo "✅ Inserted new GENOTYPE_RUN_ID: $GENOTYPE_RUN_ID" | tee -a "$LOG_FILE"
  else
    echo "✅ Found existing GENOTYPE_RUN_ID: $GENOTYPE_RUN_ID" | tee -a "$LOG_FILE"
  fi

# ---- Final Upload to Database ----
# copy Files to database

echo "Copying CSVs to database tables..." | tee -a "$LOG_FILE"

PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "\copy stock_sample FROM '$OUTPUT_STOCK_SAMPLE_CSV' WITH CSV HEADER;"
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "\copy sample_varietyset FROM '$OUTPUT_SAMPLE_VARIETYSET' WITH CSV HEADER;"
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "\copy snp_feature FROM '$OUTPUT_SNP_FEATURE' WITH CSV HEADER;"
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "\copy variant_variantset FROM '$OUTPUT_VARIANT_VARIANTSET' WITH CSV HEADER;"
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "\copy snp_featureloc FROM '$OUTPUT_SNP_FEATURELOC' WITH CSV HEADER;"

echo "✅ CSVs copied to database tables." | tee -a "$LOG_FILE"



echo "Upload completed at $(date)" >> "$LOG_FILE"


