# Workflow

<p align="center">
<img src="images/workflow.png" alt="VCF Matrix Diagram" style="width: 50%;"/>
</p>
---

## 📥 Input

- `VCF` file: `169.vcf.gz`

<img src="images/Input.png" alt="input" style="width: 70%;"/>
---

## 🔄 VCF to Matrix

### 🔧 Script: `vcftomatrix.sh`

```bash
invcf=169.vcf.gz
pref=169

mkdir -p $pref

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' $invcf | \
  tr "|" "/" | \
  tr -d "/" | \
  sed "s:chr0\\?::" | \
  awk -v OFS="\t" '{$1 = $1 + 0; print $0}' > ${pref}/mat_vcf.txt

cut -f1-4 $pref/mat_vcf.txt > $pref/pos.txt

# The following fails for some reason (grep CHR part)
bcftools view -o - -h $invcf | \
  grep CHR | \
  tr "\t" "\n" | \
  tail -n +10 > $pref/sample_list.txt
```

### 📤 Output Files

- `pos.txt`
- `sample_list.txt`
- `mat_vcf.txt`

<img src="images/Output.png" alt="output" style="width: 70%;"/>

---

## 🔄 Matrix to HDF5

### ▶️ Command

```bash
./make_HDF_dataset.sh <input-dir> <output-prefix>
```

Example:

```bash
./make_HDF_dataset.sh /home/scratch3/169 169
```

### 📤 Output File

- `169_transp.h5`

<img src="images/vmat.png" alt="hdf5" style="width: 60%;"/>

<p align="center">
<img src="images/hdf5Ui.png" alt="hdf5" style="width: 90%;"/>
</p>

---
