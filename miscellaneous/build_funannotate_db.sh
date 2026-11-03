#!/usr/bin/env bash
# build_funannotate_db.sh: build a new funannotate database snapshot for FunFlux, record the
# versions of the databases in it, and pack it with its checksum, ready to upload to Zenodo.
#
# What this is for
#   funannotate, the gene prediction and annotation program FunFlux runs, needs its own reference
#   database: UniProtKB/Swiss-Prot, Pfam, InterPro, dbCAN, MEROPS, MIBiG, the Gene Ontology, a
#   repeat-protein library, curated gene names and products (gene2product), BUSCO gene sets, and
#   the AUGUSTUS parameters of 166 trained species. `funannotate setup` builds that database by
#   downloading every part from its own server, and it reads the download addresses from
#   funannotate's GitHub repository each time it runs, so one moved or blocked server breaks it.
#   FunFlux therefore never runs the setup: rule funannotate_db downloads a frozen copy of a
#   finished database (a "snapshot") from Zenodo (links.funannotate_link in the config), checks it
#   against its checksum and unpacks it, and rule normalize_funannotate_db repairs the absolute
#   paths stored inside it. This script makes the next snapshot. You run it by hand when the
#   snapshot should be refreshed; the workflow never calls it.
#
# Usage
#   build_funannotate_db.sh OUT_DIR [BUSCO_SET ...]
#     OUT_DIR    folder for everything this script writes; created if missing. Use a folder of
#                its own, not the directories.funannotate_db of a FunFlux run.
#     BUSCO_SET  one or more of funannotate 1.8.17's BUSCO sets. The fungal ones are fungi,
#                microsporidia, dikarya, ascomycota, pezizomycotina, eurotiomycetes,
#                sordariomycetes, saccharomycetes, saccharomycetales and basidiomycota.
#                Default: dikarya, the set FunFlux trains its gene finders on unless
#                parameters.prediction.busco_db names another one.
#   Example:
#     conda activate /path/to/the/funannotate/environment
#     bash miscellaneous/build_funannotate_db.sh /data/funannotate_db_build dikarya
#
# What you need
#   * FunFlux's funannotate environment (funannotate 1.8.17), activated with `conda activate`.
#     Putting its bin/ folder on the PATH is not enough: activation also sets
#     AUGUSTUS_CONFIG_PATH, where the setup finds the trained AUGUSTUS species. Use the
#     environment Snakemake built from workflow/envs/funannotate.yaml (a folder under your
#     --conda-prefix, or under .snakemake/conda/; `snakemake --list-conda-envs` with the same
#     options as your runs prints which one), or make one from the pinned package list:
#         conda create -n funannotate_db_build --file workflow/envs/funannotate.linux-64.pin.txt
#         conda activate funannotate_db_build
#     The database should be built by the same environment that later reads it in FunFlux, so
#     that the DIAMOND and HMMER versions that write its search databases are the ones that read
#     them.
#   * Internet access to the servers in funannotate's download list (EBI, geneontology.org, the
#     dbCAN and MIBiG servers, GitHub). The build of the evening of 8 October 2026 took 10
#     minutes: about 4 for the setup, 4.5 for step 5 and 1.5 for step 6.
#   * About 9 GB of free space in OUT_DIR, on a file system that allows hard links (step 5 makes
#     a hard-linked copy inside OUT_DIR; the usual Linux file systems do, some network or USB
#     drive formats do not). The build of 8 October 2026 made a 6.5 GB database folder and a
#     2.1 GB archive (1 GB = 1,000,000,000 bytes).
#   * An OUT_DIR whose path holds none of the characters # [ * \ (step 5 uses the path as a
#     search text for sed, which reads these four as instructions; the script checks this
#     first).
#   Commands used: funannotate (with the wget, DIAMOND, HMMER and AUGUSTUS of its environment),
#   tar, gzip, sha256sum, and the standard tools of every Linux system (coreutils, awk, find,
#   grep, sed).
#
# What it writes into OUT_DIR
#   NAME below stands for funannotate_db_<funannotate version>_<YYYY-MM-DD>, for example
#   funannotate_db_1.8.17_2026-10-08; the date is the day the build started.
#     funannotate_db_<YYYY.MM.DD>/  the database as funannotate setup wrote it, with the real
#                                   paths of this computer inside. It stays after packing; you
#                                   may delete it once the archive is uploaded.
#     NAME.tar.gz        the archive to upload. Its only top-level entry is the folder above. It
#                        names no path of this computer (the stored paths read
#                        /path/to/funannotate_db) and no user or group (owner and group 0).
#     NAME.sha256        its SHA-256 checksum, one line "<hash>  NAME.tar.gz". Upload it too.
#     NAME.versions.tsv  version, date and number of records of every database in the archive.
#     NAME.setup.log     funannotate's own log of the setup. It holds local paths, so it is kept
#                        out of the archive.
#     NAME.console.log   everything the setup printed, wget's download messages included.
#   While step 5 runs, two temporary entries exist as well, and the script removes them itself,
#   also when it stops on an error: .pack_NAME/ (a hidden hard-linked copy of the database
#   folder) and NAME.tar.gz.part (the archive being written).
#
# Steps
#   1. check the arguments and the environment, and work out the names;
#   2. run funannotate setup;
#   3. check that every database arrived;
#   4. write the version table;
#   5. replace the stored paths with a placeholder in a hard-linked copy, pack the copy with
#      owner and group 0, and write the checksum;
#   6. read the archive back (checksum, top-level entry, owners) and print what to do next.
#
# Exit status
#   0  the archive and its checksum were written and checked.
#   1  stopped by one of the script's own checks; the line that starts with ERROR (or, when
#      OUT_DIR is missing, the usage line) says why.
#   Any other failure (for example tar on a full disk, or an OUT_DIR you may not write to)
#   stops the script at once (set -e, below) with that command's own message and exit status,
#   which need not be 1, and with no ERROR line.

# -e: stop at the first command that fails. -u: using a variable that was never set is an error
# (instead of silently giving empty text). -o pipefail: a pipeline (a | b) fails when any of its
# commands fails, not only the last one.
set -euo pipefail


# ───────────── 1. Arguments, environment and names ─────────────

# The funannotate version FunFlux runs (workflow/envs/funannotate.yaml). The archive name
# carries it, because the layout of the database belongs to the version that wrote it.
REQUIRED_VERSION="1.8.17"

# The ten databases that `funannotate setup -i all` installs next to the BUSCO sets. funannotate
# predict and annotate read nine of them; the tenth, busco_outgroups, serves funannotate compare.
# A complete snapshot holds all ten, each with a line in the database index,
# funannotate-db-info.txt. (A bash "array": a list of words, written in round brackets.)
REQUIRED_DATABASES=(merops uniprot dbCAN pfam repeats go mibig interpro busco_outgroups gene2product)

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 OUT_DIR [BUSCO_SET ...]   (default BUSCO set: dikarya; see the top of this script)" >&2
    exit 1
fi
out_dir=$1
shift

# Whatever follows OUT_DIR on the command line names the BUSCO sets. "${busco_sets[@]}" hands
# each set to a command as a word of its own.
if [ "$#" -gt 0 ]; then
    busco_sets=("$@")
else
    busco_sets=(dikarya)
fi
for busco_set in "${busco_sets[@]}"; do
    if [ "$busco_set" = "all" ]; then
        echo "ERROR: name the BUSCO sets one by one instead of 'all' (which downloads all 28 of them), so that step 3 can check that each one arrived." >&2
        exit 1
    fi
done

if ! command -v funannotate > /dev/null; then
    echo "ERROR: funannotate is not on the PATH. Activate FunFlux's funannotate environment first (see 'What you need' at the top of this script)." >&2
    exit 1
fi

# funannotate prints its version as "funannotate v1.8.17" and then exits with status 1 (a
# quirk of funannotate), so "|| true" keeps set -e from stopping the script here. Only the
# printed text is compared.
version_line=$(funannotate --version 2> /dev/null || true)
if [ "$version_line" != "funannotate v${REQUIRED_VERSION}" ]; then
    echo "ERROR: this script builds the database with funannotate ${REQUIRED_VERSION}, the version FunFlux runs, but 'funannotate --version' printed '${version_line}'. Activate FunFlux's funannotate environment." >&2
    exit 1
fi

# `conda activate` sets AUGUSTUS_CONFIG_PATH to the environment's AUGUSTUS config/ folder. The
# setup copies the trained species from its species/ subfolder, and crashes when the variable is
# not set. ${AUGUSTUS_CONFIG_PATH:-} reads as empty text when the variable is not set (set -u
# would stop the script otherwise).
if [ -z "${AUGUSTUS_CONFIG_PATH:-}" ] || [ ! -d "${AUGUSTUS_CONFIG_PATH}/species" ]; then
    echo "ERROR: AUGUSTUS_CONFIG_PATH is not set, or holds no species/ folder. Activate the funannotate environment with 'conda activate': putting its bin/ folder on the PATH is not enough." >&2
    exit 1
fi

# Step 2 makes funannotate download every file with wget, which the environment contains.
if ! command -v wget > /dev/null; then
    echo "ERROR: wget is not on the PATH. FunFlux's funannotate environment contains it: activate that environment." >&2
    exit 1
fi

# The day of the build, in the two spellings the names use. Every name is worked out once,
# here, so a build that runs past midnight still gives all its files the same date.
build_date=$(date +%Y-%m-%d)                     # e.g. 2026-10-08
folder_date=$(echo "$build_date" | tr '-' '.')   # e.g. 2026.10.08

# The database folder is named like the one in the December 2025 snapshot
# (funannotate_db_2025.12.15). The checksum file must be called NAME.sha256 when the archive is
# NAME.tar.gz: FunFlux works out the address of the checksum file from the archive's address
# (00_common.smk, section 10), so the two names have to stay paired.
db_folder="funannotate_db_${folder_date}"
name="funannotate_db_${REQUIRED_VERSION}_${build_date}"
archive="${name}.tar.gz"
checksum="${name}.sha256"
versions_table="${name}.versions.tsv"
setup_log="${name}.setup.log"
console_log="${name}.console.log"
# Step 5's temporary entries: a hidden hard-linked copy of the database folder, and the
# archive while tar is still writing it.
staging=".pack_${name}"
archive_part="${archive}.part"

# Work inside OUT_DIR from here on. Three things need that: the checksum line must name the
# archive without a folder in front (rule funannotate_db runs `sha256sum -c` inside the folder
# that holds both files); the archive must hold the database folder as its only top-level entry;
# and funannotate setup writes its log, funannotate-setup.log, into the folder it is started
# from, so starting it here keeps that log (it holds local paths) out of the database folder and
# out of the archive.
mkdir -p "$out_dir"
cd "$out_dir"
out_dir=$(pwd)    # the same folder as an absolute path, for funannotate and for the messages

# A folder left behind by an earlier, failed setup cannot be finished: funannotate writes the
# database index only at the very end, and a second setup on the half-built folder crashes.
if [ -e "$db_folder" ]; then
    echo "ERROR: ${out_dir}/${db_folder} already exists. If an earlier build failed, delete that folder (and any .pack_funannotate_db_* folder and *.tar.gz.part file next to it); otherwise choose another OUT_DIR. Then run this script again." >&2
    exit 1
fi
# Never overwrite an archive or a checksum file that may already have been uploaded.
if [ -e "$archive" ] || [ -e "$checksum" ]; then
    echo "ERROR: ${out_dir} already holds ${archive} or ${checksum}. Move them away, or choose another OUT_DIR." >&2
    exit 1
fi
# Leftovers of an earlier run that was killed during step 5 (a killed script cannot clean up
# after itself), from any day: a hidden .pack_funannotate_db_* copy, which holds a second name
# for every file of an old database folder (so deleting that folder alone frees no space), or
# a half-written *.tar.gz.part. find lists them (-maxdepth 1: only OUT_DIR itself).
leftovers=$(find . -maxdepth 1 \( -name '.pack_funannotate_db_*' -o -name '*.tar.gz.part' \) | sort)
if [ -n "$leftovers" ]; then
    echo "ERROR: ${out_dir} holds leftovers of an interrupted build; delete them (and the database folder of that build, if it is still there), then run this script again:" >&2
    echo "$leftovers" >&2
    exit 1
fi
# Step 5 uses the folder's path as the search text of sed, which reads # [ * and \ as
# instructions rather than as letters, so such a path would not be found there. Checked now,
# before the downloads, rather than after them. pwd -P is the same path with every symbolic
# link resolved, which step 5 searches for too.
for path_to_check in "$out_dir" "$(pwd -P)"; do
    case "$path_to_check" in
        *'#'* | *'['* | *'*'* | *'\'*)
            echo "ERROR: the path of OUT_DIR (${path_to_check}) contains one of the characters # [ * \\, which step 5 cannot use. Choose another OUT_DIR." >&2
            exit 1
            ;;
    esac
done

echo "[1/6] funannotate ${REQUIRED_VERSION}; BUSCO sets: ${busco_sets[*]}; output folder: ${out_dir}"


# ───────────── 2. Build the database with funannotate setup ─────────────
# -i all   every database funannotate knows: the ten in REQUIRED_DATABASES, and the BUSCO sets
#          named with -b. The trained AUGUSTUS species are always copied in, from
#          AUGUSTUS_CONFIG_PATH.
# -b       the BUSCO sets.
# -d       the new database folder, as an absolute path (setup creates the folder, and writes
#          this path into the database index and into every trained species' info.json).
# --wget   download with wget instead of funannotate's own downloader, which uses Python's
#          urllib: geneontology.org refuses that downloader (HTTP 403 on 8 October 2026), and
#          with this environment's Python 3.8 it does not follow the HTTP 308 redirects some
#          servers send. wget handles both.
# funannotate reads the download addresses from its GitHub repository (master branch) at run
# time, so the build fetches what those addresses serve today, including the newest UniProt,
# Pfam, InterPro and Gene Ontology releases. Everything the setup prints goes to NAME.console.log.
echo "[2/6] Running funannotate setup into ${db_folder} (about 4 minutes; its output goes to ${console_log}) ..."
if ! funannotate setup -i all -b "${busco_sets[@]}" -d "${out_dir}/${db_folder}" --wget > "$console_log" 2>&1; then
    echo "ERROR: funannotate setup failed. Its last lines:" >&2
    tail -n 20 "$console_log" >&2
    echo "Full output: ${out_dir}/${console_log}" >&2
    # Setup creates its own log and the database folder only after it has read its arguments.
    # When it stops before that (for example on a misspelt BUSCO set name), neither exists,
    # so each of these two hints is printed only when its file or folder is there.
    if [ -e funannotate-setup.log ]; then
        echo "Setup's own log: ${out_dir}/funannotate-setup.log" >&2
    fi
    if [ -e "$db_folder" ]; then
        echo "Delete ${out_dir}/${db_folder} before you run this script again." >&2
    fi
    exit 1
fi
mv funannotate-setup.log "$setup_log"


# ───────────── 3. Check that every database arrived ─────────────
# funannotate setup does not notice a failed download: it starts wget and does not look at the
# result. wget then leaves an empty or partial file behind, and the database can end up in the
# index with 0 records. So each required database must be listed in the index, with more than
# 0 records, and the file or folder the index names for it must exist.
db_index="${db_folder}/funannotate-db-info.txt"
if [ ! -s "$db_index" ]; then
    echo "ERROR: ${out_dir}/${db_index} is missing or empty: funannotate setup did not finish. See ${console_log}." >&2
    exit 1
fi
for database in "${REQUIRED_DATABASES[@]}"; do
    # The index has one tab-separated line per database: name, type, path, version, date,
    # number of records, and an md5 checksum of one downloaded file (for uniprot and pfam, the
    # short file that names the current release: reldate.txt, Pfam.version.gz), which
    # `funannotate setup -u` compares with the file on the server to spot a new release. awk
    # prints the line whose first column is this database's name (nothing, if there is no
    # such line).
    index_line=$(awk -F '\t' -v wanted="$database" '$1 == wanted' "$db_index")
    if [ -z "$index_line" ]; then
        echo "ERROR: ${database} is not listed in ${out_dir}/${db_index}. See ${console_log}." >&2
        exit 1
    fi
    db_path=$(echo "$index_line" | cut -f 3)
    db_records=$(echo "$index_line" | cut -f 6)
    if [ ! -e "$db_path" ]; then
        echo "ERROR: the index lists ${database} at ${db_path}, but that file or folder is missing. See ${console_log}." >&2
        exit 1
    fi
    # The number of records must be a whole number above 0. [[ text =~ pattern ]] compares the
    # text with a pattern; ^[0-9]+$ means "one or more digits and nothing else", so an empty or
    # garbled column stops the script too (a plain [ ... -le 0 ] test would let it through).
    if ! [[ "$db_records" =~ ^[0-9]+$ ]] || [ "$db_records" -eq 0 ]; then
        echo "ERROR: the index gives ${database} '${db_records}' records instead of a number above 0, so its download probably failed. See ${console_log}." >&2
        exit 1
    fi
done

# Each BUSCO set is unpacked into a folder named after it (FunFlux looks for the same folder
# when it checks parameters.prediction.busco_db, 00_common.smk, section 9d). dataset.cfg is the
# set's description file, which step 4 reads.
for busco_set in "${busco_sets[@]}"; do
    if [ ! -f "${db_folder}/${busco_set}/dataset.cfg" ]; then
        echo "ERROR: the BUSCO set ${busco_set} is missing or incomplete (no ${out_dir}/${db_folder}/${busco_set}/dataset.cfg). See ${console_log}." >&2
        exit 1
    fi
done

# The trained species: one folder per species under trained_species/. The folder itself is
# tested first: on a missing folder find fails, and set -e would then stop the script with
# find's own message instead of this explanation.
if [ ! -d "${db_folder}/trained_species" ]; then
    echo "ERROR: ${out_dir}/${db_folder}/trained_species is missing: funannotate setup copied no trained species. See ${console_log}." >&2
    exit 1
fi
species_count=$(find "${db_folder}/trained_species" -mindepth 1 -maxdepth 1 -type d | wc -l)
if [ "$species_count" -eq 0 ]; then
    echo "ERROR: ${out_dir}/${db_folder}/trained_species holds no species. See ${console_log}." >&2
    exit 1
fi
echo "[3/6] All ${#REQUIRED_DATABASES[@]} databases, the BUSCO sets and ${species_count} trained species are in place."


# ───────────── 4. Write the version table ─────────────
# NAME.versions.tsv: two lines starting with # that describe the build, then a table with
#   * one row per database from the index: name, type, version, date and number of records
#     (the index without its path and md5 columns). For repeats, mibig and busco_outgroups the
#     date is the day of the build, because funannotate 1.8.17 records no release date for them;
#   * one row per BUSCO set, from the set's dataset.cfg;
#   * one row for the trained species: the AUGUSTUS version they were copied from, the day of
#     the build, and how many there are.
# `augustus --version` prints "AUGUSTUS (3.5.0) is a gene prediction tool." on its first line:
# awk keeps the second word of that line and tr deletes the brackets, leaving 3.5.0.
augustus_version=$(augustus --version 2>&1 | awk 'NR == 1 {print $2}' | tr -d '()')
{
    echo "# funannotate database snapshot ${name} (folder ${db_folder}), built on ${build_date}"
    echo "# command: funannotate setup -i all -b ${busco_sets[*]} -d ${db_folder} --wget (funannotate ${REQUIRED_VERSION})"
    printf 'name\ttype\tversion\tdate\trecords\n'
    cut -f 1,2,4,5,6 "$db_index"
    for busco_set in "${busco_sets[@]}"; do
        # dataset.cfg holds lines such as name=dikarya_odb9, creation_date=2016-02-13 and
        # number_of_BUSCOs=1312; awk prints the text after the "=" of the line it is asked for.
        dataset_cfg="${db_folder}/${busco_set}/dataset.cfg"
        set_version=$(awk -F '=' '$1 == "name" {print $2}' "$dataset_cfg")
        set_date=$(awk -F '=' '$1 == "creation_date" {print $2}' "$dataset_cfg")
        set_buscos=$(awk -F '=' '$1 == "number_of_BUSCOs" {print $2}' "$dataset_cfg")
        printf 'busco_%s\tbusco_set\t%s\t%s\t%s\n' "$busco_set" "$set_version" "$set_date" "$set_buscos"
    done
    printf 'trained_species\taugustus\t%s\t%s\t%s\n' "$augustus_version" "$build_date" "$species_count"
} > "$versions_table"
echo "[4/6] Version table written: ${versions_table}"


# ───────────── 5. Remove the traces of this computer, pack the folder, write its checksum ─────────────
# The archive is published, so it must not reveal the computer that built it. Two things would:
#   * funannotate setup writes the absolute path of the database folder (the -d of step 2) into
#     the database index (funannotate-db-info.txt, third column) and into the info.json of every
#     trained species (the "path" of its AUGUSTUS parameters). In the archive these paths become
#     the placeholder /path/to/funannotate_db. Nothing is lost: on another computer the real
#     paths would be wrong anyway, so whoever unpacks the archive has to rewrite them. FunFlux's
#     rule normalize_funannotate_db does that for any stored path (it keeps the part from
#     trained_species/ on, or the end part that exists in the new folder), and so does the
#     manual recipe in docs/getting-started/databases.md ("Using the snapshot without
#     FunFlux"), which reads the old prefix from the index.
#   * tar records in every entry the name and number of the user and group who own the file.
#     --owner=0 --group=0 --numeric-owner records 0 and no name instead. Whoever unpacks the
#     archive as an ordinary user owns the unpacked files anyway.
# The database folder itself keeps its real paths, so that it still works as a database where
# it is (docs/getting-started/databases.md, "Building a new snapshot"). The placeholder goes
# into a copy made of hard links (cp -al): a second name for every file, which takes no extra
# disk space. sed -i writes each file it edits anew under the same name, so in the copy the
# edited files become new files, and the originals in the database folder stay as they were.
placeholder="/path/to/funannotate_db"
echo "[5/6] Replacing the stored paths with ${placeholder} in a hard-linked copy of ${db_folder} ..."
# From here to the end of this step, the two temporary entries (named in step 1) are removed
# whenever the script stops, on an error too: trap runs the command in quotes when the script
# exits. Removing the copy removes only the second names; the files of the database folder
# stay. ${staging:?} stops rm instead of running it if the name were ever empty. The trap is
# switched off again ("trap - EXIT") once the archive has its final name.
trap 'rm -rf "${staging:?}" "${archive_part:?}"' EXIT
mkdir "$staging"
cp -al "$db_folder" "${staging}/${db_folder}"

# The files that hold the paths: the index and the info.json of every trained species. find
# lists those; "read -r" takes its output one line (one file) at a time into the bash list.
path_files=("${staging}/${db_index}")
while IFS= read -r info_json; do
    path_files+=("$info_json")
done < <(find "${staging}/${db_folder}/trained_species" -mindepth 2 -maxdepth 2 -name info.json)

# The folder's path in two spellings: as handed to setup (out_dir, from pwd), and with every
# symbolic link resolved (pwd -P), in case setup stored that one. sed replaces both, the longer
# one first: one spelling can end in the other (a link /scratch pointing to /gpfs/scratch), and
# replacing the shorter one first would then leave the start of the longer one in the files.
# (Step 1 refused the characters sed would misread; a "." is a wildcard for sed but also matches
# itself, so it changes nothing here.) ${#name} is the length of the text in name.
physical_out_dir=$(pwd -P)
logical_prefix="${out_dir}/${db_folder}"
physical_prefix="${physical_out_dir}/${db_folder}"
if [ "${#physical_prefix}" -ge "${#logical_prefix}" ]; then
    build_prefixes=("$physical_prefix" "$logical_prefix")
else
    build_prefixes=("$logical_prefix" "$physical_prefix")
fi
for build_prefix in "${build_prefixes[@]}"; do
    sed -i "s#${build_prefix}#${placeholder}#g" "${path_files[@]}"
done

# No file of the copy may still name the database folder's real path, in either spelling: the
# exact text setup wrote and sed replaced. (Not OUT_DIR alone: a short OUT_DIR such as /data
# also occurs in the downloaded data itself, in README files of trained species.) grep -r
# searches every file (all 6.5 GB, about a minute), -l lists the files with a match, -F reads the
# paths as plain text. grep's exit status says: 0 = found (stop), 1 = nothing found (go on),
# other = grep failed.
grep_status=0
leftover=$(grep -rlF -e "${out_dir}/${db_folder}" -e "${physical_out_dir}/${db_folder}" "${staging}/${db_folder}") || grep_status=$?
if [ "$grep_status" -eq 0 ]; then
    echo "ERROR: these files still name ${out_dir}/${db_folder}; the archive was not written. Look at them in ${out_dir}:" >&2
    # The paths grep printed start with the copy's folder, which the trap removes; sed cuts that
    # part off, so the list names the same files in the database folder, which stays.
    echo "$leftover" | sed "s#^${staging}/##" >&2
    exit 1
elif [ "$grep_status" -ne 1 ]; then
    echo "ERROR: grep could not search ${out_dir}/${staging} (exit status ${grep_status}); the archive was not written." >&2
    exit 1
fi

# tar -c creates an archive, -z compresses it with gzip, -f names the archive file, and -C makes
# tar read from the copy, so the folder (under its own name) is the archive's only top-level
# entry. Rule funannotate_db unpacks with `tar -xzf ... --strip-components=1`, which drops
# exactly that one level, so the files land directly in directories.funannotate_db.
# tar writes NAME.tar.gz.part first, and the file gets its final name only when tar has finished,
# so an interrupted run never leaves a half-written archive under the name you would upload.
echo "[5/6] Packing ${db_folder} into ${archive} (about 4 minutes) ..."
tar --owner=0 --group=0 --numeric-owner -czf "$archive_part" -C "$staging" "$db_folder"
mv "$archive_part" "$archive"

# Done with the copy: remove it now and switch the trap off.
rm -rf "${staging:?}"
trap - EXIT

# The checksum line reads "<hash>  NAME.tar.gz", with no folder in front of the file name,
# because the command runs in the folder that holds the archive.
sha256sum "$archive" > "$checksum"


# ───────────── 6. Read the archive back, and say what comes next ─────────────
# Does the archive match its checksum, and is the database folder its only top-level entry?
# tar -t lists every path in the archive; cut keeps the part before the first "/", and sort -u
# keeps each different value once.
echo "[6/6] Checking the archive ..."
if ! sha256sum -c "$checksum"; then
    echo "ERROR: ${archive} does not match ${checksum}." >&2
    exit 1
fi
top_level=$(tar -tzf "$archive" | cut -d '/' -f 1 | sort -u)
if [ "$top_level" != "$db_folder" ]; then
    echo "ERROR: the archive's top-level entries are not just ${db_folder}:" >&2
    echo "$top_level" >&2
    exit 1
fi
echo "Top-level entry of ${archive}: ${top_level}/"

# Every entry must record owner and group 0 (step 5). With --numeric-owner, the second column
# of `tar -tv` is "owner number/group number" without names; awk prints that column.
owners=$(tar --numeric-owner -tvzf "$archive" | awk '{print $2}' | sort -u)
if [ "$owners" != "0/0" ]; then
    echo "ERROR: the archive's entries record other owners than 0/0:" >&2
    echo "$owners" >&2
    exit 1
fi
echo "Owner and group of every entry: 0/0; stored paths: ${placeholder}"

# Sizes, both in decimal units (1 GB = 1,000,000,000 bytes) like the notes at the top: the
# folder as the disk space it takes (du -s -B 1 counts it in bytes), the archive as its exact
# number of bytes. numfmt turns a number of bytes into a short size such as 2.1G. It rounds up
# unless told otherwise (2,140,223,992 bytes would print as 2.2G), hence --round=nearest.
# (du -h would count in powers of 1024 and print 6.1G for the 6.5 GB folder.)
folder_bytes=$(du -s -B 1 "$db_folder" | cut -f 1)
folder_size=$(numfmt --to=si --round=nearest "$folder_bytes")
archive_bytes=$(stat -c %s "$archive")
archive_size=$(numfmt --to=si --round=nearest "$archive_bytes")

cat << EOF

Snapshot ready in ${out_dir}:
  ${archive}    ${archive_size} (${archive_bytes} bytes)
  ${checksum}
  ${versions_table}
  ${setup_log}, ${console_log}    for your records (local paths: do not upload)
  ${db_folder}/    ${folder_size} on disk (real local paths inside; not needed once the archive is uploaded)

Next steps:
  1. Upload these files to Zenodo as a NEW VERSION of the record that links.funannotate_link
     points at, so that the record's concept DOI stays the same (the new version gets a record
     number of its own):
         ${archive}
         ${checksum}
         ${versions_table}    (optional)
  2. In the FunFlux config, set
         links.funannotate_link: https://zenodo.org/records/<new record number>/files/${archive}?download=1
     FunFlux finds ${checksum} in the same record by its name.
  3. Point directories.funannotate_db at a NEW, empty folder. Never reuse a folder that holds an
     older snapshot: FunFlux may keep using the old snapshot there, or unpack the new one over it
     and leave old files behind, and other workflows may still read that folder.
EOF
