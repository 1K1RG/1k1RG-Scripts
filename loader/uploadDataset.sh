#!/bin/bash





# ---- Arguments to file Inputs ------

DB="$1"
VARIANTSET_NAME="$2"
SAMPLE_FILE="$3"
POS_FILE="$4"
ORGANISM_NAME="$5"
ORGANISM_GENUS="$5"
ORGANISM_SPECIES="$5"
H5_file="$6"


OUTPUT_DIR="$VARIANTSET_NAME"
mkdir -p "$OUTPUT_DIR"

timenow=$(date '+%Y-%m-%d %H:%M:%S')


# ---- Usage Check ----
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <sample_file> <pos_file>"
  exit 1
fi


if [[ ! -f "$SAMPLE_FILE" ]]; then
  echo "❌ Error: File '$SAMPLE_FILE' not found."
  exit 1
fi

if [[ ! -f "$POS_FILE" ]]; then
  echo "❌ Error: File '$POS_FILE' not found."
  exit 1
fi

# ---- Configuration ----
DB_NAME="20k"
DB_USER="iricadmin"
DB_HOST="localhost"
DB_PORT="5432"
DB_PASSWORD="user12345"
CSV_DIR="/path/to/csv/files"
LOG_FILE="upload_log.txt"

# ---- Variable Definition ----
TERM_NAME='whole plant'
CV_NAME='plant_anatomy'



CVTERM_SNP='SNP'
CVTERM_CHR='chromosome'
VARIANT_TYPE_CV='sequence'  # or use 'SO' or your CV name

OUTPUT_STOCK_SAMPLE_CSV="$OUTPUT_DIR/StockSample.csv"
OUTPUT_SAMPLE_VARIETYSET="$OUTPUT_DIR/SampleVarietySet.csv"
OUTPUT_SNP_FEATURE="$OUTPUT_DIR/SnpFeature.csv"
OUTPUT_VARIANT_VARIANTSET="$OUTPUT_DIR/Variant_variantset.csv"
OUTPUT_SNP_FEATURELOC="$OUTPUT_DIR/Snp_FeatureLoc.csv"

# ---- Start Upload ----
echo "Starting bulk upload at $timenow" > "$LOG_FILE"

# -- INITIALIZE MAX 
MAX_DB_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -t -A -q -c \
	  "SELECT COALESCE(MAX(db_id), 0) FROM db;")

  
  MAX_DB_ID=$(echo "$MAX_DB_ID" | xargs)
  MAX_DB_ID=$((MAX_DB_ID + 1))

  PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -c \
	  "SELECT setval('db_db_id_seq', $MAX_DB_ID, false);"

MAX_DBXREF_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -t -A -q -c \
	  "SELECT COALESCE(MAX(dbxref_id), 0) FROM dbxref;")

  
  MAX_DBXREF_ID=$(echo "$MAX_DBXREF_ID" | xargs)
  MAX_DBXREF_ID=$((MAX_DBXREF_ID + 1))

  PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -c \
	  "SELECT setval('dbxref_dbxref_id_seq', $MAX_DBXREF_ID, false);"


# ---- CVTERM Lookup or Insert ----
echo "Looking up CVTERM__PLANT_ID for '$TERM_NAME' from '$CV_NAME'..." | tee -a "$LOG_FILE"

CVTERM__PLANT_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
SELECT cvterm_id FROM cvterm JOIN cv ON cvterm.cv_id = cv.cv_id 
WHERE cvterm.name = '$TERM_NAME' AND cv.name = '$CV_NAME' LIMIT 1;
")

CVTERM__PLANT_ID=$(echo $CVTERM__PLANT_ID | xargs)

if [[ -z "$CVTERM__PLANT_ID" ]]; then
  echo "⚠️ CVTERM__PLANT_ID not found, inserting cv and cvterm..." | tee -a "$LOG_FILE"

  # Ensure CV exists
  CV_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
  INSERT INTO cv(name) 
  VALUES ('$CV_NAME') 
  ON CONFLICT (name) DO UPDATE SET name=EXCLUDED.name 
  RETURNING cv_id;
  ")

  CV_ID=$(echo $CV_ID | xargs)

  # Insert cvterm
  CVTERM__PLANT_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
  INSERT INTO cvterm(name, cv_id, is_obsolete) 
  VALUES ('$TERM_NAME', $CV_ID, false)
  RETURNING CVTERM__PLANT_ID;
  ")
fi

CVTERM__PLANT_ID=$(echo $CVTERM__PLANT_ID | xargs)
echo "✅ GOT CVTERM__PLANT_ID: $CVTERM__PLANT_ID" | tee -a "$LOG_FILE"

# ---- Organism Lookup or Insert ----
echo "Looking up organism_id for '$ORGANISM_NAME'..." | tee -a "$LOG_FILE"

ORGANISM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
SELECT organism_id FROM organism WHERE common_name = '$ORGANISM_NAME' LIMIT 1;
")

if [[ -z "$ORGANISM_ID" ]]; then
  echo "⚠️ organism_id not found, inserting organism..." | tee -a "$LOG_FILE"

  ORGANISM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
  INSERT INTO organism (genus, species, common_name)
  VALUES ('$ORGANISM_GENUS', '$ORGANISM_SPECIES', '$ORGANISM_NAME')
  RETURNING organism_id;
  ")
fi

ORGANISM_ID=$(echo $ORGANISM_ID | xargs)
echo "✅ GOT ORGANISM_ID: $ORGANISM_ID" | tee -a "$LOG_FILE"

# FEATURE TERM

FEATURE_CVTERM_CHR_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
	SELECT cvterm.cvterm_id 
	FROM cvterm 
	JOIN cv ON cv.cv_id = cvterm.cv_id 
	WHERE cv.name = '$VARIANT_TYPE_CV' AND cvterm.name = '$CVTERM_CHR'
	LIMIT 1;
	")

FEATURE_CVTERM_CHR_ID=$(echo $FEATURE_CVTERM_CHR_ID | xargs)

# ---- DB Lookup or Insert ----
echo "Looking up db_id for '$DB'..." | tee -a "$LOG_FILE"


DB_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
SELECT db_id FROM db WHERE name = '$DB' LIMIT 1;
")

if [[ -z "$DB_ID" ]]; then
  echo "⚠️ db_id not found, inserting db..." | tee -a "$LOG_FILE"

  DB_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
  INSERT INTO db(name)
  VALUES ('$DB')
  RETURNING db_id;")

   DB_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
    SELECT db_id FROM db WHERE name = '$DB' LIMIT 1;
  ")
fi

DB_ID=$(echo $DB_ID | xargs)
echo "✅ GOT DB_ID: $DB_ID" | tee -a "$LOG_FILE"

# ---- VariantSet Lookup or Insert ----

echo "Looking up variantset_id for variantset name '$VARIANTSET_NAME'..." | tee -a "$LOG_FILE"

VARIANTSET_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
SELECT variantset_id FROM variantset WHERE name = '$VARIANTSET_NAME' LIMIT 1;")

if [[ -z "$VARIANTSET_ID" ]]; then
  echo "⚠️ variantset_id not found." | tee -a "$LOG_FILE"
  echo "Looking up CVTERM_PLANT_ID for variant type '$CVTERM_SNP' in CV '$VARIANT_TYPE_CV'..." | tee -a "$LOG_FILE"

	VARIANTTYPE_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
	SELECT cvterm.cvterm_id
	FROM cvterm 
	JOIN cv ON cv.cv_id = cvterm.cv_id 
	WHERE cv.name = '$VARIANT_TYPE_CV' AND cvterm.name = '$CVTERM_SNP'
	LIMIT 1;
	")

	VARIANTTYPE_ID=$(echo $VARIANTTYPE_ID | xargs)

	if [[ -z "$VARIANTTYPE_ID" ]]; then
	  echo "❌ Could not find CVTERM__PLANT_ID for variant type '$CVTERM_SNP' in '$VARIANT_TYPE_CV'" | tee -a "$LOG_FILE"
	  exit 1
	fi

   echo "✅ GOT VARIANTTYPE_ID: $VARIANTTYPE_ID" | tee -a "$LOG_FILE"

   echo "Adding variantset_id." | tee -a "$LOG_FILE"
  
   echo "Fetching current max variant_set_id..." | tee -a "$LOG_FILE"

  MAX_VARIANT_SET_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -t -A -q -c \
	  "SELECT COALESCE(MAX(variantset_id), 0) FROM variantset;")

	MAX_VARIANT_SET_ID=$(echo "$MAX_VARIANT_SET_ID" | xargs)
	VARIANT_SET_ID=$((MAX_VARIANT_SET_ID + 1))

	PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -c \
	  "SELECT setval('variantset_variantset_id_seq', $VARIANT_SET_ID, false);"


  VARIANTSET_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
  INSERT INTO variantset(name, variant_type_id, organism_id)
  VALUES ('$VARIANTSET_NAME', $VARIANTTYPE_ID, $ORGANISM_ID)
  RETURNING variantset_id;
  ")
fi

VARIANTSET_ID=$(echo $VARIANTSET_ID | xargs)
echo "✅ GOT VARIANTSET_ID: $VARIANTSET_ID" | tee -a "$LOG_FILE"

# ---- Platform Lookup or Insert ----

echo "Looking up platform_id for '$PLATFORM_NAME'..." | tee -a "$LOG_FILE"

PLATFORM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
SELECT platform_id FROM platform WHERE variantset_id=$VARIANTSET_ID and db_id=$DB_ID LIMIT 1;
")

if [[ -z "$PLATFORM_ID" ]]; then
  echo "⚠️ platform_id not found, inserting platform..." | tee -a "$LOG_FILE"

  MAX_PLATFORM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -q -c "SELECT COALESCE(MAX(platform_id), 0) FROM platform;" | xargs)

	# Update the sequence (replace sequence name with actual name)
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('platform_platform_id_seq1', $MAX_PLATFORM_ID + 1, false);"



  PLATFORM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
  INSERT INTO platform (variantset_id, db_id, genotyping_method_id)
  VALUES ($VARIANTSET_ID, $DB_ID, NULL)
  RETURNING platform_id;
  ")

fi

PLATFORM_ID=$(echo $PLATFORM_ID | xargs)
echo "✅ GOT PLATFORM_ID: $PLATFORM_ID" | tee -a "$LOG_FILE"

#---- ADD GenotypeRun

GENOTYPE_RUN_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
SELECT genotype_run_id FROM genotype_run WHERE platform_id=$PLATFORM_ID and data_location=$H5_file LIMIT 1;
")

if [[ -z "$GENOTYPE_RUN_ID" ]]; then
  echo "⚠️ GENOTYPE_RUN_ID not found, inserting GENOTYPE_RUN_ID..." | tee -a "$LOG_FILE"

  MAX_GENOTYPE_RUN_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -q -c "SELECT COALESCE(MAX(genotype_run_id), 0) FROM genotype_run;" | xargs)

	# Update the sequence (replace sequence name with actual name)
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('genotype_run_genotype_run_id_seq', $MAX_GENOTYPE_RUN_ID + 1, false);"


  GENOTYPE_RUN_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
  INSERT INTO genotype_run (platform_id, date_performed, data_location)
  VALUES ($PLATFORM_ID, $(date) , $H5_file)
  RETURNING genotype_run_id;
  ")

fi

PLATFORM_ID=$(echo $PLATFORM_ID | xargs)
echo "✅ GOT PLATFORM_ID: $PLATFORM_ID" | tee -a "$LOG_FILE"


# ---- Add Sample to database ----

echo "Starting sample loading..." | tee "$LOG_FILE"

MAX_STOCK_SAMPLE_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -q -c "SELECT COALESCE(MAX(stock_sample_id), 0) FROM stock_sample;" | xargs)

echo "🔢 Max stock_sample_id: $MAX_STOCK_SAMPLE_ID"

	# Update the sequence (replace sequence name with actual name)
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('stock_sample_stock_sample_id_seq', $MAX_STOCK_SAMPLE_ID + 1, false);"

STOCK_SAMPLE_ID=$((MAX_STOCK_SAMPLE_ID +1))

echo "✅ (STOCK_SAMPLE) Sequence updated to start at $((STOCK_SAMPLE_ID))"


echo "stock_sample_id, stock_id,dbxref_id,hdf5_index,tmp_oldstock_id" > "$OUTPUT_STOCK_SAMPLE_CSV"

# --- Get max sample_varietyset_id ---
MAX_SAMPLE_VARIETYSET_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -q -c "SELECT COALESCE(MAX(sample_varietyset_id), 0) FROM sample_varietyset;" | xargs)

SAMPLE_VARIETYSET_ID=$((MAX_SAMPLE_VARIETYSET_ID +1))


# --- Update the sequence ---
PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('sample_run_sample_run_id_seq', $((SAMPLE_VARIETYSET_ID)), false);
"

# --- Get max sample_varietyset_id ---
MAX_STOCK_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" \
  -d "$DB_NAME" -t -A -q -c "SELECT COALESCE(MAX(stock_id), 0) FROM stock;" | xargs)

MAX_STOCK_ID=$((MAX_STOCK_ID +1))

PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -c "
  SELECT setval('stock_stock_id_seq', $((MAX_STOCK_ID)), false);
"


echo "🔢 Max sample_varietyset_id: $SAMPLE_VARIETYSET_ID"


echo "✅ (SAMPLE VARIETY_SET)  Sequence updated to start at $((SAMPLE_VARIETYSET_ID))"


echo "sample_varietyset_id, stock_sample_id, db_id,hdf5_index" > "$OUTPUT_SAMPLE_VARIETYSET"


COUNTER=0
HDF_COUNTER=0

while IFS= read -r LINE || [ -n "$LINE" ]; do
  echo "Processing line: $LINE" | tee -a "$LOG_FILE"

     # ---- Add / Check for Stock Record ----
  LINE=$(printf "%s" "$LINE" | sed "s/'/''/g" | sed 's/\\/\\\\/g') 	
		
  STOCK_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
    SELECT stock_id FROM stock
    WHERE uniquename = '$LINE'
	  AND name = '$LINE'
      AND type_id = $CVTERM__PLANT_ID
      AND organism_id = $ORGANISM_ID
    LIMIT 1;
  " | xargs)

  if [[ -z "$STOCK_ID" ]]; then
    echo "⚠️ Stock '$LINE' not found. Inserting..." | tee -a "$LOG_FILE"

    STOCK_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
      INSERT INTO stock (name, uniquename, type_id, organism_id)
      VALUES ('$LINE', '$LINE', $CVTERM__PLANT_ID, $ORGANISM_ID)
      RETURNING stock_id")

    echo "✅ Inserted new stock_id: $STOCK_ID at $HDF_COUNTER" | tee -a "$LOG_FILE"
  else
    echo "✅ Found existing stock_id: $STOCK_ID at $HDF_COUNTER" | tee -a "$LOG_FILE"
  fi

  
  
  # --- DBXREF INSERT OR SELECT ---
	echo "Looking up dbxref_id for db_id=$DB_ID and accession='$LINE'..." | tee -a "$LOG_FILE"

	DBXREF_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
	  SELECT dbxref_id FROM dbxref
	  WHERE db_id = $DB_ID AND accession = '$LINE'
	  LIMIT 1;
	" | xargs)

	if [[ -z "$DBXREF_ID" ]]; then
	  echo "⚠️ dbxref not found. Inserting..." | tee -a "$LOG_FILE"

	  DBXREF_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
		INSERT INTO dbxref (db_id, accession, version)
		VALUES ($DB_ID, '$LINE', 1)
		RETURNING dbxref_id;
	  " | xargs)

	  echo "✅ Inserted new dbxref_id: $DBXREF_ID" | tee -a "$LOG_FILE"
	else
	  echo "✅ Found existing dbxref_id: $DBXREF_ID" | tee -a "$LOG_FILE"
	fi
	
  
  echo "$STOCK_SAMPLE_ID, $STOCK_ID,$DBXREF_ID, $HDF_COUNTER, $HDF_COUNTER" >> "$OUTPUT_STOCK_SAMPLE_CSV"
  
  echo "$SAMPLE_VARIETYSET_ID, $STOCK_SAMPLE_ID,$DB_ID, $HDF_COUNTER" >> "$OUTPUT_SAMPLE_VARIETYSET"
  
  
  ((STOCK_SAMPLE_ID++))
   ((SAMPLE_VARIETYSET_ID++))
  ((HDF_COUNTER++))
 

done < "$SAMPLE_FILE"

FEATURE_CVTERM_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
	SELECT cvterm_id 
	FROM cvterm 
	JOIN cv ON cv.cv_id = cvterm.cv_id 
	WHERE cv.name = '$VARIANT_TYPE_CV' AND cvterm.name = '$CVTERM_CHR'
	LIMIT 1;
	")
	
	FEATURE_CVTERM_ID=$(echo "$FEATURE_CVTERM_ID" | xargs)

MAX_SNPFEATURE_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -t -A -q -c \
	  "SELECT COALESCE(MAX(snp_feature_id), 0) FROM snp_feature;")

MAX_SNPFEATURE_ID=$(echo "$MAX_SNPFEATURE_ID" | xargs)
SNP_FEATURE_ID=$((MAX_SNPFEATURE_ID + 1))

PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -c \
	  "SELECT setval('snp_feature_snp_feature_id_seq', $SNP_FEATURE_ID, false);"





MAX_VARIANT_VARIANTSET_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -t -A -q -c \
	  "SELECT COALESCE(MAX(variant_variantset_id), 0) FROM variant_variantset;")

MAX_VARIANT_VARIANTSET_ID=$(echo "$MAX_VARIANT_VARIANTSET_ID" | xargs)
VARIANT_VARIANTSET_ID=$((MAX_VARIANT_VARIANTSET_ID + 1))

PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -c \
	  "SELECT setval('variant_variantset_variant_variantset_id_seq', $VARIANT_VARIANTSET_ID, false);"


MAX_FEATURELOC_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -t -A -q -c \
	  "SELECT COALESCE(MAX(snp_featureloc_id), 0) FROM snp_featureloc;")

MAX_FEATURELOC_ID=$(echo "$MAX_FEATURELOC_ID" | xargs)
SNP_FEATURELOC_ID=$((MAX_FEATURELOC_ID + 1))

PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" \
	  -U "$DB_USER" -d "$DB_NAME" -c \
	  "SELECT setval('snp_featureloc_snp_featureloc_id_seq', $SNP_FEATURELOC_ID, false);"


counter=0;

echo "snp_feature_id, variantset_id" > "$OUTPUT_SNP_FEATURE"

echo "variant_variantset_id, variant_feature_id, variant_type_id, variantset_id, hdf5_index" > "$OUTPUT_VARIANT_VARIANTSET"

echo "snp_featureloc_id, snp_feature_id, organism_id, srcfeature_id, snp_feature_id, position, refcall" > "$OUTPUT_SNP_FEATURELOC"

echo "Starting position loading..." | tee -a "$LOG_FILE"
{
	while IFS=$'\t' read -r chrom posit refc rest; do
		if [[ "$chrom" != "$POS" ]]; then
			echo "reading $POS_FILE"
			echo "chromosome: $chrom"
			NAME="chr${chrom,,}"
			if [[ "$chrom" =~ ^[0-9]+$ ]]; then
			    UNIQUENAME=$(printf "chr%02d" "$chrom")
			else
			    UNIQUENAME="chr${chrom,,}"
			fi

			
			echo "Inserting new feature for $NAME $UNIQUENAME $counter..."
			
			echo "NAME: $NAME"
			echo "UNIQUENAME: $UNIQUENAME"
			echo "FEATURE_CVTERM_ID: $FEATURE_CVTERM_ID"
			echo "ORGANISM_ID: $ORGANISM_ID"

			
			FEATURE_ID=$(PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -t -A -q -c "
				SELECT feature_id 
				FROM feature
				WHERE name ilike '$NAME' 
					AND uniquename ilike '$UNIQUENAME'
					AND type_id = $FEATURE_CVTERM_ID
					AND organism_id = $ORGANISM_ID
				LIMIT 1;
				")

			FEATURE_ID=$(echo "$FEATURE_ID" | xargs)
			echo "FEATURE_ID is..$FEATURE_ID"
			POS=$chrom
		fi
		
		echo "$SNP_FEATURE_ID, $VARIANTSET_ID" >&3
		
		echo "$VARIANT_VARIANTSET_ID,$SNP_FEATURE_ID,\\N,$VARIANTSET_ID,$counter" >&4
		
		FINAL_POS=$((posit - 1))
		
		echo "$SNP_FEATURELOC_ID,$SNP_FEATURE_ID,$ORGANISM_ID,$FEATURE_ID,$FINAL_POS,$refc" >&5
		
		((SNP_FEATURE_ID++))
		((VARIANT_VARIANTSET_ID++))
		((SNP_FEATURELOC_ID++))
		#echo "adding $counter"
		

		((counter++))
	done < "$POS_FILE"
} 3>>"$OUTPUT_SNP_FEATURE" 4>>"$OUTPUT_VARIANT_VARIANTSET" 5>>"$OUTPUT_SNP_FEATURELOC"

echo "Loading CSVs into database..." | tee -a "$LOG_FILE"

PGPASSWORD=$DB_PASSWORD psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" <<EOF
\copy stock_sample FROM '$OUTPUT_STOCK_SAMPLE_CSV' WITH (FORMAT csv, HEADER true)
\copy sample_varietyset FROM '$OUTPUT_SAMPLE_VARIETYSET' WITH (FORMAT csv, HEADER true)
\copy snp_feature FROM '$OUTPUT_SNP_FEATURE' WITH (FORMAT csv, HEADER true)
\copy variant_variantset FROM '$OUTPUT_VARIANT_VARIANTSET' WITH (FORMAT csv, HEADER true, NULL '\N');
\copy snp_featureloc FROM '$OUTPUT_SNP_FEATURELOC' WITH (FORMAT csv, HEADER true, NULL '\N')
EOF

echo "CSV import completed at $timenow" | tee -a "$LOG_FILE"


echo "Upload completed at $timenow" >> "$LOG_FILE"
