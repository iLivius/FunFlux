# Reference databases

FunFlux needs reference databases in addition to its conda environments. You provide three of
them before the first run: NCBI `core_nt` with its taxonomy files (not read in `contigs` mode
without the screen), the eggNOG-mapper data and InterProScan. The workflow, BUSCO and Medaka fetch
the rest, on the computer that runs Snakemake, except BUSCO with `auto`, which downloads inside
each job. `snakemake download_all` fetches all of it in advance, except BUSCO's data with `auto`,
as described in
[Preparation, on a computer with internet](../reference/running.md#cluster-preparation).

The table lists every database with its configuration key, when that key is required, who
provides the data, its version and its size.

| Database | Config key | Key required | Provided by | Version | Size |
|---|---|---|---|---|---|
| [NCBI `core_nt` and taxonomy files](#core-nt) | `directories.blast_db`; `parameters.nt_version` names the database | with the contamination and mitochondrion screen: `illumina`, `longread`, `hybrid`, and `contigs` with `parameters.contigs.screen: true` | you | NCBI's current build; tested with the build of 11 July 2026 | 264 GB |
| [eggNOG-mapper data](#eggnog) | `directories.eggnog_db` | yes | you | eggNOG 5.0.2 | 48 GB |
| InterProScan member databases | `files.iprscan` | yes | you, with [InterProScan](external-tools.md#interproscan) | 5.78-109.0 | 7.0 GB download, 50 GB installed |
| [funannotate database](#funannotate-db) | `directories.funannotate_db`, `links.funannotate_link` | yes | the workflow, when the folder holds no database | `funannotate_db_1.8.17_2026-10-08` | 2.1 GB archive; 6.1 GB unpacked, 8.1 GB with the archive |
| [UNITE](#unite) | `links.unite_its_link`, or `directories.unite_db` | the link, unless `directories.unite_db` is set | the workflow, once per output folder; or your copy | UNITE release for eukaryotes, 19.02.2025 | 35.7 MB download, 207 MB unpacked |
| [antiSMASH databases](#antismash) | `directories.antismash_db` | no | the workflow, once per output folder; or your copy | as downloaded by antiSMASH 8.0.4 | 9.4 GB |
| [Tiberius model](#tiberius-model) | `links.tiberius_model_link`, `parameters.tiberius.model_sha256` | the link, unless `parameters.tiberius.mode` is `off` | the workflow, once per output folder | the `fungi` model of 8 February 2026 | 30 MB download, 31 MB unpacked |
| [Telomere repeat database](#telomere-db) | `links.telomere_db_link`, `parameters.contig_summary.telomere_db_sha256` | yes | the workflow, once per output folder | `telo_database_v2.tsv` of tidk, at the commit the template link names | 0.6 MB |
| [ChroQueTas and the FungAMR database](#chroquetas) | `links.chroquetas_link`, `parameters.amr.chroquetas_script_sha256`, `parameters.amr.chroquetas_db_sha256` | only with `parameters.amr.chroquetas: on` (off by default) | the workflow, once per output folder | ChroQueTas 1.0.1 and its database `db-v20250811` (57 species) | 29 kB script, 137 kB archive; 2.5 MB unpacked |
| [BUSCO lineage datasets](#busco) | `parameters.busco.lineage` | no (default `auto`) | BUSCO, in every BUSCO job with `auto`; the workflow, once per output folder, for a named lineage | odb12 | per archive, e.g. eukaryota_odb12 29 MB, hypocreales_odb12 275 MB |
| [PhiX genome](#phix) | `links.phix_link` | in `illumina` and `hybrid` | the workflow, once per output folder (`download_phix`) | GCF_000819615.1 | 2 kB |
| [Medaka model](#medaka) | `parameters.longread.medaka_model`, `parameters.hybrid.medaka_model` | no (default `auto`) | Medaka, during the model check, the first time a model is used | named, or worked out from the reads | not measured |

Sizes on disk are as `du -h` prints them, in powers of 1024. Download and archive sizes are the
byte counts the servers list, in decimal units. The databases you keep yourself take about 380 GB
with a local antiSMASH copy, and 264 GB of that is `core_nt`, which only the screen reads.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **taxonomy files** | the files `taxdb.btd` and `taxdb.bti`, which give BLAST its organism names, and `nodes.dmp` and `names.dmp`, the NCBI taxonomy BlobTools reads; all beside `core_nt` |
    | **NOTE** | a line starting `[FunFlux] NOTE` that reports something to check without stopping the run; read them after every dry run |
    | **funannotate database** | the frozen snapshot of funannotate's data, such as Swiss-Prot, Pfam, the BUSCO set and stored AUGUSTUS models, downloaded once from Zenodo into `directories.funannotate_db` |
    | **stored AUGUSTUS model** | the AUGUSTUS parameters trained earlier for one species, kept in `trained_species/` of the funannotate database; named in column 4 of the annotation parameter table, they are used as they are |
    | **path repair** | rule `normalize_funannotate_db`: it rewrites the absolute paths stored in the funannotate database to your folder, in place and without a lock, so start one run first |
    | **launch folder** | the folder Snakemake is started from, normally the FunFlux repository; relative paths are read from it, and `.snakemake/` is created there |
    | **BUSCO set** | the BUSCO genes funannotate looks for to train SNAP, GlimmerHMM and, without a stored AUGUSTUS model, AUGUSTUS: `dikarya` (`dikarya_odb9`) unless `parameters.prediction.busco_db` names another |
    | **lineage dataset** | the BUSCO gene set of one taxonomic group, such as `hypocreales_odb12`, for the completeness check; with `lineage: auto`, BUSCO chooses one by placing the genome |

## NCBI `core_nt` and taxonomy files { #core-nt }

The contamination and mitochondrion screen searches every contig against this database with
`blastn -task megablast` (rule `blast_contigs`). BlobTools reads the NCBI taxonomy from the same
folder (rule `blob_json`), so `directories.blast_db` must hold three kinds of files:

| Files | From | Used for |
|---|---|---|
| `core_nt.nal`, `core_nt.00.*`, `core_nt.01.*`, ... | `update_blastdb.pl core_nt` | the search; `parameters.nt_version` (default `core_nt`) is the database name |
| `taxdb.btd`, `taxdb.bti` | `taxdb.tar.gz` | organism names in the BLAST table; without them the name columns read `N/A` |
| `nodes.dmp`, `names.dmp` | `taxdump.tar.gz` | the taxonomy tree and names BlobTools uses to assign a taxon to each contig |

`parameters.nt_version: nt` selects the full `nt` database instead, if that is what the folder
holds. These commands download `core_nt` and the taxonomy files:

```bash
conda create -n blast -c conda-forge -c bioconda blast=2.16.0   # provides update_blastdb.pl
conda activate blast
mkdir -p /path/to/blast_db && cd /path/to/blast_db
update_blastdb.pl --decompress --num_threads 4 core_nt
wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar -zxf taxdb.tar.gz
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && tar -zxf taxdump.tar.gz
```

`update_blastdb.pl` checks every archive against the `.md5` file NCBI publishes beside it, and
`--decompress` unpacks each archive and deletes it. The build of 11 July 2026 has 80 volumes
(252 GB) and 12 GB of database-wide index files, and the taxonomy files add 0.8 GB. FunFlux v2
does not read `nucl_gb.accession2taxid` (13 GB), which the v1 instructions also downloaded.

Two commands check the database and the taxonomy files:

```bash
BLASTDB=/path/to/blast_db blastdbcmd -db core_nt -info | head -4
ls /path/to/blast_db/taxdb.btd /path/to/blast_db/taxdb.bti /path/to/blast_db/nodes.dmp /path/to/blast_db/names.dmp
```

For the build of 11 July 2026, the first command prints:

```text
Database: Core nucleotide BLAST database
	128,460,233 sequences; 947,384,893,387 total bases

Date: Jul 11, 2026  1:17 AM	Longest sequence: 98,200,793 bases
```

!!! warning "A dry run that passes does not show that the volumes and `taxdb` are in place"

    A missing folder is noticed before any job starts, with a NOTE. A missing `nodes.dmp` or
    `names.dmp` is noticed in a dry run with the screen on, which stops with
    `MissingInputException in rule blob_json` because both files are inputs of `blob_json`.
    Missing volumes or a missing `taxdb` are noticed only when BLAST runs, which first reads them.

## eggNOG-mapper data { #eggnog }

FunFlux runs eggNOG-mapper 2.1.13, which reads the data of eggNOG 5.0.2. The data for
eggnog-mapper 3.x (eggNOG 7) do not work with it, so download the data with eggNOG-mapper
2.1.13 itself:

```bash
conda create -n eggnog-mapper -c conda-forge -c bioconda eggnog-mapper=2.1.13
conda activate eggnog-mapper
mkdir -p /path/to/eggnog_db
download_eggnog_data.py --data_dir /path/to/eggnog_db -y
```

The script fetches `eggnog.db`, `eggnog.taxa.db` and the DIAMOND database `eggnog_proteins.dmnd`
from `eggnogdb.embl.de` with `wget`, which the environment contains, and `-y` answers its
questions with yes. The Pfam, MMseqs2 and HMMER databases are downloaded only on request (`-P`,
`-M`, `-H`). FunFlux does not use them, because it runs `emapper.py -m diamond`.

To check the download, list the folder:

```bash
ls -lh /path/to/eggnog_db
# eggnog.db 39G, eggnog_proteins.dmnd 8.7G, eggnog.taxa.db 266M, eggnog.taxa.db.traverse.pkl 6.4M
```

FunFlux prints a NOTE before any job starts when the folder does not exist. The files themselves
are first read when `eggnog_annotation` runs, after gene prediction.

## funannotate database snapshot { #funannotate-db }

funannotate needs its own database. FunFlux does not build one with `funannotate setup` but takes
a frozen snapshot from Zenodo, `funannotate_db_1.8.17_2026-10-08`, record 23265189
([10.5281/zenodo.23265189](https://doi.org/10.5281/zenodo.23265189)). It was written on
8 October 2026 by `funannotate setup -i all -b dikarya --wget` with funannotate 1.8.17, the
version FunFlux runs. All versions of the funannotate database on Zenodo share the DOI
[10.5281/zenodo.18271294](https://doi.org/10.5281/zenodo.18271294), which leads to the newest.

The database holds the components below, with the versions to cite for the annotation databases.
Most of these versions are those funannotate wrote into the database index. The `dikarya_odb9`
row comes from the description file of that BUSCO set, and the row of the stored AUGUSTUS models
from the AUGUSTUS version that built the database. The Zenodo record also holds the versions as a
table, `funannotate_db_1.8.17_2026-10-08.versions.tsv`.

| Component | Version | Records |
|---|---|---|
| UniProtKB/Swiss-Prot | 2026_03 | 575,748 |
| Pfam | 38.2 | 30,134 |
| InterPro | 110.0 | 54,922 |
| Gene Ontology | 2026-07-26 | 48,340 |
| gene2product | 1.99 | 34,505 |
| dbCAN | HMMdb V11 (funannotate's index calls it 14.0) | 699 |
| MEROPS | 12.5 (the label funannotate 1.8.17 gives it) | 5,098 |
| MIBiG | 1.4 | 31,023 |
| repeat proteins | 1.0 | 11,950 |
| BUSCO outgroups | 1.0 | 8 in the index; `outgroups/` holds the sets of 6 species |
| BUSCO `dikarya_odb9` | 2016-02-13 | 1,312 BUSCOs |
| stored AUGUSTUS models, in `trained_species/` (AUGUSTUS and other gene-finder parameters) | from AUGUSTUS 3.5.0 | 166 |

??? note "Five databases are newer than in the previous snapshot (FunFlux v1)"

    The previous snapshot, `funannotate_db_1.8.15_2025-12-15` (record 18271295, named for
    funannotate 1.8.15), stays on Zenodo, and FunFlux v1 downloads it. A folder that already holds
    it keeps it. To use the new snapshot, point `directories.funannotate_db` at a new, empty
    folder, for the reasons given under [Newer snapshots](#build-snapshot).

    | Component | Previous snapshot (2025-12-15) |
    |---|---|
    | UniProtKB/Swiss-Prot | 2025_04, 573,661 records |
    | Pfam | 38.0, 25,545 records |
    | InterPro | 107.0, 49,674 records |
    | Gene Ontology | 2025-10-10, 48,196 records |
    | gene2product | 1.97, 34,471 records |
    | MEROPS | the same file, labelled 12.0 |
    | dbCAN, MIBiG, repeat proteins, BUSCO outgroups, BUSCO `dikarya_odb9`, stored AUGUSTUS models | the same |

    The other components hold the same data. Because five databases are newer, the gene names,
    products, domains and GO terms of a new annotation can differ from those made with the
    previous snapshot. The gene models can differ too, because funannotate predict also aligns
    Swiss-Prot to the genome as protein evidence.

### What the workflow does

Two steps prepare the funannotate database. It is downloaded once per database folder, and the
path repair checks its stored paths in every new output folder.

| Step (name in the log) | When | What it does |
|---|---|---|
| `funannotate_db` | only while `<funannotate_db>/funannotate-db-info.txt` is missing | downloads the archive named by `links.funannotate_link` and its `.sha256` file from the same Zenodo record, checks the archive with `sha256sum -c` and unpacks it into the folder; archive and checksum file stay there. Log: `logs/funannotate_db.log` |
| `normalize_funannotate_db` | in every new output folder | the funannotate database stores absolute paths in `funannotate-db-info.txt` and `trained_species/*/info.json`; the step rewrites those that do not point into your folder and leaves correct files alone. Summary: `00.resources/funannotate_db_paths.normalized.json` |

A checksum mismatch stops the run with
`[FunFlux] funannotate_db: ... does not match its checksum file ..., so it was NOT unpacked`.
Delete both files, the archive and its checksum file, and run again.

!!! warning "One run at a time until the paths are repaired"

    `normalize_funannotate_db` rewrites files inside the database folder, in place and without a
    lock. Two runs that start at the same time on a database whose paths are not yet repaired can
    damage it. Start one run and let it pass this step before you start the others on the same
    folder. Once the paths are repaired, the step only reads the database. Keep the folder
    writable until then.

    Reach the folder under the same real path in every run. A run that sees it under another path,
    as a second user or a second machine may, rewrites the stored paths. `funannotate-db-info.txt`
    then changes, and the next run on any other output folder repeats the repair and every
    annotation step after it.

### Downloading it by hand

For a machine without internet access at run time, these commands take the same steps as rule
`funannotate_db`:

```bash
mkdir -p /path/to/funannotate_db && cd /path/to/funannotate_db
curl -L --fail -o funannotate_db_1.8.17_2026-10-08.tar.gz \
    'https://zenodo.org/records/23265189/files/funannotate_db_1.8.17_2026-10-08.tar.gz?download=1'
curl -L --fail -o funannotate_db_1.8.17_2026-10-08.sha256 \
    'https://zenodo.org/records/23265189/files/funannotate_db_1.8.17_2026-10-08.sha256?download=1'
sha256sum -c funannotate_db_1.8.17_2026-10-08.sha256   # funannotate_db_1.8.17_2026-10-08.tar.gz: OK
tar -xzf funannotate_db_1.8.17_2026-10-08.tar.gz --strip-components=1
ls funannotate-db-info.txt dikarya trained_species
```

### Using the snapshot without FunFlux { #without-funflux }

The funannotate database also serves funannotate 1.8.17 on its own, once its stored paths are
rewritten the way `normalize_funannotate_db` does it:

```bash
mkdir -p /path/to/funannotate_db
tar -xzf funannotate_db_<version>_<date>.tar.gz --strip-components=1 -C /path/to/funannotate_db
cd /path/to/funannotate_db
OLD=$(head -n 1 funannotate-db-info.txt | cut -f 3 | xargs dirname)
sed -i "s#${OLD}#$(pwd)#g" funannotate-db-info.txt trained_species/*/info.json
export FUNANNOTATE_DB=/path/to/funannotate_db
funannotate database
```

The sed line reads the stored prefix from the database index and replaces it with the current
folder, in the index and in the `info.json` of every stored AUGUSTUS model. The folder's path
must not contain `&` or `#`, which sed would read as instructions. The last command should list
ten databases.

### Newer snapshots, and building your own { #build-snapshot }

Which snapshot a run uses depends on what the folder `directories.funannotate_db` holds:

| The folder | What FunFlux does |
|---|---|
| is empty, or has no `funannotate-db-info.txt` | downloads the snapshot `links.funannotate_link` names |
| holds a funannotate database | uses it, whatever `links.funannotate_link` names |
| holds a funannotate database that an earlier run from the same launch folder downloaded, and `links.funannotate_link` has changed since | unpacks the new archive over the old files and leaves behind those the new one does not replace ([What makes Snakemake re-run a rule](../reference/running.md#what-makes-snakemake-re-run-a-rule)) |

!!! warning "For a new snapshot, point `directories.funannotate_db` at a new, empty folder"

    As the table shows, a folder that already holds a funannotate database either keeps it or has
    the new archive unpacked over its old files. Other runs and workflows may also still read the
    old folder.

To build a database with the current UniProt, Pfam, InterPro and Gene Ontology releases yourself,
run `miscellaneous/build_funannotate_db.sh` in FunFlux's funannotate environment:

```bash
conda activate /path/to/the/funannotate/environment    # the folder --list-conda-envs prints
bash miscellaneous/build_funannotate_db.sh /path/to/funannotate_db_build dikarya
```

- Activate the environment with `conda activate`. With only its `bin/` folder on the `PATH`, the
  setup misses `AUGUSTUS_CONFIG_PATH`, which activation sets. The script stops with a message
  when that variable is not set, or when funannotate is missing or not version 1.8.17.
- The script needs network access and about 9 GB of free space in the folder named as its first
  argument, and takes about 10 minutes. It refuses to overwrite a database folder, archive or
  checksum file that already exists.
- It runs `funannotate setup -i all -b <BUSCO sets> -d <folder> --wget` and checks that the ten
  databases, the BUSCO sets and the stored AUGUSTUS models all arrived. It writes the database
  folder `funannotate_db_<YYYY.MM.DD>/` and an archive of it,
  `funannotate_db_<funannotate version>_<date>.tar.gz`, with its `.sha256` file. It also writes a
  table of the database versions (`.versions.tsv`) and two logs.
- Point `directories.funannotate_db` at the database folder it wrote. That folder holds
  `funannotate-db-info.txt`, so nothing is downloaded.
- The script ends by printing the steps to publish the archive on Zenodo. For your own runs you
  need none of them: keep the database folder and point `directories.funannotate_db` at it.

### Other BUSCO sets for gene prediction { #busco-sets }

funannotate predict trains SNAP and GlimmerHMM on the conserved genes of the BUSCO set named by
`parameters.prediction.busco_db`. It trains AUGUSTUS on them too, for the samples without a model
in column 4 of the
[annotation parameter table](../reference/configuration.md#annotation-parameter-table).

| Set | BUSCO genes | In the funannotate database | Tested with FunFlux v2 |
|---|---|---|---|
| `dikarya`, the default | 1,312 | yes, the only one | yes |
| `fungi` | 290 | no | no: its download has not been tested |

Keep `dikarya` unless the fungus lies outside the Dikarya. AUGUSTUS needs at least 200 genes to
train on, so with `fungi` a sample is safer with a stored AUGUSTUS model in column 4. When a
funannotate database is present and the named set is not, every Snakemake command stops before any
job starts and prints the command that adds the set,
`funannotate setup -i busco -b <set> -d <folder> -w`. Run it in FunFlux's funannotate
environment, activated as below.

For a fungus outside the Dikarya (Ascomycota and Basidiomycota), add a set such as `fungi` with
funannotate 1.8.17, for example from the environment FunFlux built. Run the first command while
`busco_db` is still `dikarya`, with the `--conda-prefix` of your runs, or without it if you use
none:

```bash
snakemake --sdm conda --list-conda-envs --conda-prefix /path/to/funflux_envs \
    --configfile config/config_custom.yaml          # prints the folder of funannotate.yaml
conda activate /path/to/that/folder
funannotate setup -i busco -b fungi -d /path/to/funannotate_db -w
```

If the folder does not exist yet, build the environments first, as described in
[Building the environments in advance](installation.md#build-in-advance). Then set
`parameters.prediction.busco_db: fungi`.

- Activate the environment with `conda activate`. With only its `bin/` folder on the `PATH`, setup
  stops with `KeyError: 'AUGUSTUS_CONFIG_PATH'`.
- `-i busco` limits the downloads to BUSCO sets. Setup still copies into `trained_species/` the
  AUGUSTUS species files the database lacks (it already holds all 166), rewrites
  `funannotate-db-info.txt` from the lines it can read, and writes `funannotate-setup.log` into the
  current folder.
- `-w` makes funannotate download with `wget` instead of its own downloader, which has failed on
  some of funannotate's download links: an HTTP 308 redirect on the `dikarya` link, and HTTP 403
  from geneontology.org.

Adding a set this way has been tried only with the small `microsporidia` set. The `fungi`
download has not been tested, and FunFlux v2 was tested with `dikarya` only.

## UNITE { #unite }

FunFlux classifies the ITS1 and ITS2 sequences with `vsearch --sintax` against UNITE. The default
`links.unite_its_link` is the UNITE USEARCH/UTAX/SINTAX release for eukaryotes, version 19.02.2025
(DOI 10.15156/BIO/3301246), which covers all eukaryotes, not Fungi only. FunFlux either downloads
this file or links your local copy:

| | Default: download | Local copy: `directories.unite_db` |
|---|---|---|
| Rule | `unite_db` | `unite_db_local` |
| What it does | streams the file through `gzip -dc` into `00.resources/unite/unite_its_sintax.fasta`, once per output folder | links the file into `00.resources/unite/` |
| `00.resources/unite/PROVENANCE.txt` records | the link, the download date, the md5 and the size | the folder, the date, the md5 and the size |
| Check | no checksum is verified | the run stops before any job starts when the folder holds no uncompressed `unite_its_sintax.fasta` |

A local copy serves any number of output folders. To make one, download and unpack the file
yourself:

```bash
mkdir -p /path/to/unite && cd /path/to/unite
curl -L --fail 'https://s3.hpc.ut.ee/plutof-public/original/338a1413-6039-4e00-b5cf-410346a1e366.gz' \
    | gzip -dc > unite_its_sintax.fasta
grep -c '^>' unite_its_sintax.fasta   # 266589
md5sum unite_its_sintax.fasta         # 4cacdc7f2372e0d1e863d7f9715a5843
```

```yaml
directories:
  unite_db: /path/to/unite
```

When `links.unite_its_link` is also set, a NOTE says that it is not used. Another UNITE release
in USEARCH/UTAX/SINTAX format can replace the default through either key, and its classifications
can then differ from those of the default release.

Downloaded or linked, the file is then filtered. Rule `unite_identified` keeps the records whose
taxonomy names a genus in `00.resources/unite/unite_its_sintax.identified.fasta`, and SINTAX reads
that file unless `parameters.its.unite_identified_only` is `false`. See
[The UNITE reference filter](../analysis/its-taxonomy.md#the-unite-reference-filter).

## antiSMASH databases { #antismash }

The antiSMASH databases end up in `00.resources/antismash/databases/` of the output folder, either
downloaded there or linked from a copy you hold:

| | Default: download | Local copy: `directories.antismash_db` |
|---|---|---|
| Rule | `antismash_db` | `antismash_db_local` |
| What it does | runs `download-antismash-databases` of antiSMASH 8.0.4 into `00.resources/antismash/databases/`, once per output folder, from `dl.secondarymetabolites.org` and `ftp.ebi.ac.uk` (Pfam 35.0); the folder, 9.4 GB, stays for later runs in the same output folder | fills `00.resources/antismash/databases/` with one symbolic link per entry of your folder; nothing is copied, and FunFlux does not write into your folder |

A local copy saves 9.4 GB in every output folder and serves any number of them. To make one:

```bash
conda create -n antismash -c conda-forge -c bioconda antismash=8.0.4
conda activate antismash
download-antismash-databases --database-dir /path/to/antismash_db
ls /path/to/antismash_db
# as-js  clusterblast  clustercompare  comparippson  knownclusterblast  mite  nrps_pks  pfam  resfam  tigrfam
```

The downloader ends with `Pre-building all databases...` and `done.`. Use a database made by
antiSMASH 8.0.4, the version FunFlux runs.

```yaml
directories:
  antismash_db: /path/to/antismash_db
```

The run stops before any job starts when the folder does not exist or has no `clusterblast/`.
Otherwise FunFlux prints
`Using the local antiSMASH database at '/path/to/antismash_db'. Nothing will be downloaded.`

## BUSCO lineage datasets { #busco }

There is nothing to download beforehand. Where BUSCO 6.0.0 gets its lineage data from
`busco-data.ezlab.org` depends on the setting:

| `parameters.busco.lineage` | BUSCO option | Downloaded | Where and when |
|---|---|---|---|
| `auto` (default) | `--auto-lineage-euk` | eukaryota_odb12 (29 MB), the eukaryote placement files (8 MB) and the lineage BUSCO places the genome in | by every BUSCO job, into `02.assembly/{sample}/eval/busco/busco_tmp`, which Snakemake deletes when the job ends. Samples never share a download, because BUSCO 6.0.0 has no lock on its download folder and silently updates a shared one |
| a lineage name, e.g. `hypocreales_odb12` | `-l hypocreales_odb12 --offline` | that dataset | once per output folder, by rule `busco_data`, into `00.resources/busco/`, on the computer that runs Snakemake; every BUSCO job reads it there |

!!! note "With `auto`, every BUSCO job needs network access"

    A named lineage needs network access once, and the downloaded dataset stays at the version in
    `00.resources/busco/` for later runs into the same output folder. BUSCO fetches the newest
    version of a dataset, so with `auto` completeness can also shift when a new version is
    published.

`busco --list-datasets` in a BUSCO 6.0.0 environment lists the names, each with its version after
a dot, such as `hypocreales_odb12.2`. Write the name without the dot and the number after it:
`hypocreales_odb12`.

The BUSCO set that funannotate trains its gene finders on (`parameters.prediction.busco_db`,
`dikarya_odb9`) is a different dataset. It is part of the funannotate database, as described
under [Other BUSCO sets for gene prediction](#busco-sets).

## Tiberius model { #tiberius-model }

The gene finder [Tiberius](../analysis/gene-prediction.md#tiberius), on by default, needs the
`fungi` model its authors trained, a set of network weights. Rule `tiberius_model` downloads it
from `links.tiberius_model_link`, on the Tiberius authors' server, once per output folder. It
checks the archive against `parameters.tiberius.model_sha256` with `sha256sum -c` and unpacks it
into `00.resources/tiberius/model/` (`model_config.json`, `model_layers.json`, `weights.h5`).
FunFlux does not ship or mirror the model, because no licence is published with it.

The address carries no version number. When the authors replace the model there, the checksum no
longer matches, the run stops before anything is unpacked, and the downloaded file is kept as
`00.resources/tiberius/model.tar.gz`. Either set `parameters.tiberius.mode: off`, or check the new
model and accept it:

```bash
sha256sum /path/to/output_dir/00.resources/tiberius/model.tar.gz
# put the 64 characters it prints into parameters.tiberius.model_sha256, then run again
```

A new checksum re-runs Tiberius for every sample, and the gene prediction after it.

## Telomere repeat database { #telomere-db }

The contig summary chooses the telomere repeat unit of every assembly from the data and checks the
units it weighs against the database that comes with tidk (Brown et al. 2025). That database,
`telo_database_v2.tsv`, is a tab-separated table of 4,051 rows, one per examined species, and a
few species appear twice. Each row gives the telomere repeat unit reported for the species, its
taxonomy from genus to superkingdom, the assembly the unit was checked on, and how well it is
supported:

| Support | Rows | Used by FunFlux |
|---|---|---|
| `strong` | 2,292 | yes |
| `weak` | 698 | yes |
| `unchecked` | 629 | no |
| `none`, with no unit | 432 | no |

FunFlux takes the units with `strong` or `weak` support, 256 distinct units, as candidate telomere
units in every run. Most rows are animals (2,996), while 903 are plants and 126 fungi. The
candidates table of each sample says how many species per kingdom carry each unit, as described on
[Contig summary](../analysis/contig-summary.md#telomere-unit).

Rule `telomere_db` downloads the file with `curl` from `links.telomere_db_link`, once per output
folder, and checks it against `parameters.contig_summary.telomere_db_sha256` with `sha256sum -c`.
It keeps the file as `00.resources/telomeres/telo_database_v2.tsv`, with a `PROVENANCE.txt` beside
it that records the source URL, the download date, the checksum and the row count. The template
address names one commit of the repository `tolkit/a-telomeric-repeat-database` on GitHub, so the
file behind it does not change. A file that does not match the checksum stops the run with
`[FunFlux] telomere_db: ...`, and the download is kept with the suffix `.failed`. To move to a
newer version of the database, change the link and the checksum together:

```bash
curl -L --fail -o telo_database_v2.tsv '<the new address>'
sha256sum telo_database_v2.tsv   # put the 64 characters into parameters.contig_summary.telomere_db_sha256
```

A new checksum re-runs `telomere_db`, then `telomere_units` and the contig summary of every
sample. There is no local-copy key: the file is 0.6 MB, and a run needs network access to
`raw.githubusercontent.com` once per output folder. The repository publishes no licence file, so
FunFlux neither ships the database nor mirrors it. It is downloaded on your machine at run time,
as the Tiberius model is. Cite it through the tidk paper and the address, as listed on
[Citation and references](../about/citation.md).

## ChroQueTas and the FungAMR database { #chroquetas }

The AMR screen, also called the resistance-mutation screen, is off by default. It runs ChroQueTas
against the FungAMR database, which lists point mutations with the evidence for or against
resistance to antifungal drugs and fungicides, as described on
[Antifungal resistance mutations](../analysis/amr.md). ChroQueTas is one bash script, and its
authors keep the database prepared for it in the same GitHub repository, `nmquijada/ChroQueTas`.
Both are licensed CC BY-NC-ND 4.0, so FunFlux does not ship or mirror them, and does not change
them.

With `parameters.amr.chroquetas: on`, rule `chroquetas_db` downloads three files with `curl` from
`links.chroquetas_link`, once per output folder, into `00.resources/chroquetas/`:

| File at the link | Kept as | Checked against |
|---|---|---|
| `bin/ChroQueTas.sh` (29,202 bytes) | `ChroQueTas.sh` | `parameters.amr.chroquetas_script_sha256` |
| `FungAMR_db/db-v20250811.tgz` (137,391 bytes) | unpacked into `FungAMR_db/`, one folder per species; the archive is then deleted | `parameters.amr.chroquetas_db_sha256` |
| `LICENSE` | `LICENSE` | its checksum is recorded in `PROVENANCE.txt` |

The rule checks the script and the archive with `sha256sum -c` before anything runs or is
unpacked. It then writes `species.txt`, which lists the 57 species folders one per line, and
`PROVENANCE.txt`, which records the three addresses, the download date, the checksums, the
species count and the licence. The template address names one commit, the one tagged `v1.0.1`,
so the files behind it do not change. A file that does not match its checksum stops the run with
`[FunFlux] chroquetas_db: ...` and is kept with the suffix `.failed`.

To move to a newer version, read the new script and its licence, then change the link and both
checksums together. FunFlux reads the files ChroQueTas 1.0.1 writes, so a version that writes them
differently stops `amr_chroquetas` with an error and needs a newer FunFlux.

```bash
BASE=https://raw.githubusercontent.com/nmquijada/ChroQueTas/<commit>
curl -L --fail -o ChroQueTas.sh "$BASE/bin/ChroQueTas.sh"
curl -L --fail -o db.tgz "$BASE/FungAMR_db/<archive name>"
sha256sum ChroQueTas.sh db.tgz   # into chroquetas_script_sha256 and chroquetas_db_sha256
```

The archive name, `db-v20250811.tgz`, is fixed in FunFlux, so a newer database under another name
needs a newer FunFlux. A new link or checksum re-runs `chroquetas_db`, then `amr_chroquetas` of
every sample. There is no local-copy key: the files are small, and a run with the AMR screen on
needs network access to `raw.githubusercontent.com` once per output folder.

## PhiX genome { #phix }

Only `illumina` and `hybrid` mode use the PhiX genome. `links.phix_link` points at NCBI's PhiX174
genome (GCF_000819615.1, 2 kB compressed), the control that Illumina libraries are spiked with.
`download_phix` fetches it with `curl` once per output folder into `00.resources/phix/`, and
`build_phix` indexes it there for bowtie2. `map_phix` then removes the read pairs of every sample
that map to it. In these modes an empty value stops the run before any job starts.

## Medaka model { #medaka }

A Medaka model is needed in `longread` and `hybrid` only, for ONT samples, unless
`medaka_model: false`. Medaka 2.2.2 ships the weights of four models. For any other model it
downloads them the first time the model is used, so the run needs network access once.

The download happens in `check_medaka_model`, before the assembler starts and on the computer
that runs Snakemake. A named model is checked and its weights fetched. With `medaka_model: auto`,
the model is worked out from the basecaller tag in the read headers, and its weights are fetched.
Medaka keeps them in its conda environment, or in `~/.medaka/data` when it cannot write there, so
polishing needs no network. See [The Medaka model](../methods/medaka-model.md).

## Pointing a key at a local copy { #local-copies }

Five database keys name a folder on your own disk. The table gives what each folder holds, what
FunFlux checks before any job starts, and what the run does with the folder:

| Key | The folder holds | Checked before any job starts | What the run does |
|---|---|---|---|
| `directories.blast_db` | the `core_nt` files, `taxdb.btd`, `taxdb.bti`, `nodes.dmp`, `names.dmp` | NOTE when the folder is missing; with the screen on, the dry run stops without `nodes.dmp` or `names.dmp` | reads it |
| `directories.eggnog_db` | the files `download_eggnog_data.py` writes | NOTE when the folder is missing | reads it |
| `directories.funannotate_db` | the unpacked funannotate database, or nothing yet | the `parameters.prediction.busco_db` set, when a funannotate database is present | downloads the database if it is missing; rewrites the paths inside it |
| `directories.antismash_db` | `clusterblast/` and the other database folders | the run stops when the folder is missing or has no `clusterblast/` | links its entries into `00.resources/antismash/databases/` |
| `directories.unite_db` | `unite_its_sintax.fasta`, uncompressed | the run stops when the file is missing | links the file into `00.resources/unite/` |

FunFlux turns each of these paths into an absolute one, measured from the launch folder. Write
them as absolute paths. Every key is listed on [Configuration](../reference/configuration.md).
