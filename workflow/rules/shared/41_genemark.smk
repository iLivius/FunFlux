# FunFlux v2.1.0 — GeneMark-ES licence check (rules/shared/41_genemark.smk)
#
# GeneMark-ES is an optional gene finder that funannotate can use during gene prediction. It is
# licensed per user: the licence is a key file, ~/.gm_key, that EXPIRES (GeneMark's download
# page says the licence is valid for one year from the date of download). The trouble is WHEN
# the expiry shows. gmes_petap.pl, the program funannotate calls, only checks that the key file
# EXISTS; the expiry is found out later, inside GeneMark's own compiled programs, when
# funannotate predict finally reaches the GeneMark step, which is well into the
# prediction. The message is "Your license period has ended", and the whole prediction has been
# wasted by then. A key that worked a few months earlier can be dead by the next
# run, so the date of the key file does not tell when the key expires.
#
# This module adds one small rule that makes the failure appear within a minute instead. It runs
# GeneMark-ES, with the same options funannotate uses, on a 200 kb slice of the sample's own genome
# (GeneMark-ES finishes such a slice in about a minute when the licence is valid) and stops the run
# with a clear message if the licence or the installation is broken. funannotate_prediction takes the
# result as an input, so the gene prediction, and everything after it, cannot start before this
# check has passed.
#
# Where it sits: in the annotation stage, between repeat masking and gene prediction, once per
# sample ("{sample}" in the file names is a Snakemake wildcard: the rule runs once for every
# sample, with the sample name filled in).
#
#   contigs_mask.fasta (CONTIGS_MASK, from repeat_masking) -> genemark_check
#       -> 04.annotation/{sample}/genemark_check.txt (GENEMARK_CHECK) -> funannotate_prediction
#
# Only defined when GeneMark is switched on (GENEMARK_ENABLED in 00_common.smk, which is true when
# directories.genemark_dir is set and parameters.genemark.mode is not off; 00_common.smk also
# checks then, before any job starts, that gmes_petap.pl and ~/.gm_key exist). With GeneMark off
# the rule does not exist, and funannotate_prediction waits for nothing. That input of
# funannotate_prediction comes from genemark_gate in 00_common.smk, a function Snakemake calls for
# each job; with GeneMark off it returns an empty list.
#
# Names used here, all defined in 00_common.smk: GENEMARK_ENABLED, CONTIGS_MASK, GENEMARK_CHECK,
# GENEMARKDIR (directories.genemark_dir), GENEMARK_KEY (the path of ~/.gm_key, used only in the
# message) and LOGS.
#
# conda: the path is written relative to THIS file (workflow/rules/shared/), so
# "../../envs/funannotate.yaml" lands on workflow/envs/funannotate.yaml.


if GENEMARK_ENABLED:

    # ── genemark_check — does GeneMark-ES run, with a valid licence? ──
    # Takes in:    CONTIGS_MASK from rule repeat_masking: the repeat-masked genome that funannotate
    #              predict will receive.
    # Does:        copies the first 200,000 bytes of it (usually the start of scaffold_1, the
    #              longest contig) into a temporary folder, which is deleted when the job ends.
    #              The job works inside that folder, because GeneMark writes its working files
    #              into the current folder. It starts gmes_petap.pl exactly the way funannotate
    #              predict does (see "Which Perl" below), on the slice, on one core, with the
    #              options funannotate predict gives GeneMark at funannotate's default settings,
    #              which FunFlux does not change:
    #                --ES              GeneMark's self-training mode: it learns the gene model from
    #                                  the sequence itself
    #                --fungus          adds the model of intron branch points, meant for fungi
    #                --max_intron 3000 the longest intron allowed, in bases
    #                --soft_mask 2000  lowercase letters are repeats; lowercase stretches longer
    #                                  than 2000 bases are left out of the gene search
    #              Then it reads GeneMark's output:
    #                * Output that mentions the licence (the word "license", in any case), or that
    #                  GeneMark could not start ("can't locate", Perl's message for a missing
    #                  module; "bad interpreter", the system's message for a Perl named in the
    #                  first line that does not exist; "command not found"; "no such file"), is a
    #                  failure: the rule stops with a [FunFlux] message, on the screen and in the
    #                  log, that says what to do.
    #                * Any other problem on the slice (for example, too little sequence to train on in
    #                  a very small genome) is NOT a licence problem and is only noted in the log:
    #                  the real run decides whether the genome suits GeneMark.
    #              The log holds GeneMark's exit status, its whole output and the outcome.
    # Produces:    GENEMARK_CHECK, a small text file with the date (UTC) and the outcome. Its only
    #              job is to exist: it is the proof funannotate_prediction waits for.
    # Consumed by: funannotate_prediction, through genemark_gate in 00_common.smk.
    #
    # The Conda environment is funannotate's, the environment funannotate predict runs GeneMark
    # in. PATH is cut down to that environment plus the system tools for the same reason as in
    # funannotate_prediction: so that no GeneMark installation that happens to sit elsewhere on
    # the PATH is picked up by accident. This rule calls gmes_petap.pl by its full path, and
    # gmes_petap.pl finds GeneMark's other programs in its own folder. GENEMARK_PATH is
    # therefore not set here: it is the variable funannotate predict reads to find GeneMark
    # (funannotate_prediction sets it), and no file of the GeneMark-ES folder reads it.
    #
    # Which Perl. funannotate predict runs gmes_petap.pl as a program, not as `perl
    # gmes_petap.pl`, so the first line of the script (the "#!" line) chooses the Perl. GeneMark
    # ships its scripts with "#!/usr/bin/perl", the computer's own Perl, which often lacks the
    # modules GeneMark needs (YAML, Hash::Merge, MCE, Parallel::ForkManager, ...); funannotate
    # predict then stops late, at the GeneMark step, with "Can't locate ... in @INC". This check
    # starts the script the same way, so the same Perl answers and the same failure shows
    # within a minute. The fix, which the message names, is GeneMark's own script
    # change_path_in_perl_scripts.pl: run once in the GeneMark folder as
    #   perl change_path_in_perl_scripts.pl "/usr/bin/env perl"
    # it makes every GeneMark script start with "#!/usr/bin/env perl", which picks the first
    # Perl on the PATH: in FunFlux, the Perl of the funannotate environment, which holds those
    # modules. Installing the modules into the computer's Perl works too.
    #
    # v1 had no such check (new in v2).
    rule genemark_check:
        input:
            contigs = CONTIGS_MASK,
        output:
            check = GENEMARK_CHECK,
        params:
            genemark_dir = GENEMARKDIR,
            key_file = GENEMARK_KEY,
        conda:
            "../../envs/funannotate.yaml"
        log:
            LOGS + "/genemark_check_{sample}.log"
        shell:
            """
            work_dir=$(mktemp -d)
            trap 'rm -rf "$work_dir"' EXIT
            cd "$work_dir"

            head -c 200000 {input.contigs} > test_slice.fasta

            export PATH="$CONDA_PREFIX/bin:/usr/bin:/bin"

            # The exit status is read below instead of being allowed to stop the rule at once,
            # because a failure for a harmless reason (a slice too small to train on) must not
            # stop the run.
            status=0
            "{params.genemark_dir}/gmes_petap.pl" \
                --ES \
                --fungus \
                --max_intron 3000 \
                --soft_mask 2000 \
                --cores 1 \
                --sequence test_slice.fasta > gmes_output.txt 2>&1 || status=$?

            echo "GeneMark-ES test run on a 200 kb slice, exit status $status. Output follows." > {log}
            cat gmes_output.txt >> {log}

            if grep -qi "license" gmes_output.txt; then
                echo "[FunFlux] GeneMark-ES refused to run because of its licence (key file {params.key_file}). The key has probably expired: request a new one from the GeneMark site, or set parameters.genemark.mode: off. See {log}." | tee -a {log} >&2
                exit 1
            fi
            if grep -qiE "can't locate|bad interpreter" gmes_output.txt; then
                echo "[FunFlux] GeneMark-ES could not start from {params.genemark_dir}: the Perl named in the first line of gmes_petap.pl (GeneMark ships it as /usr/bin/perl, the computer's own) is missing or lacks a Perl module GeneMark needs. funannotate starts GeneMark the same way, so its gene prediction would fail too. Fix it once, in that folder: perl change_path_in_perl_scripts.pl '/usr/bin/env perl' (GeneMark's scripts then use the Perl of FunFlux's funannotate environment, which has the modules), or install the missing modules into that Perl, or set parameters.genemark.mode: off. See {log}." | tee -a {log} >&2
                exit 1
            fi
            if grep -qiE "command not found|no such file" gmes_output.txt; then
                echo "[FunFlux] GeneMark-ES could not start from {params.genemark_dir} (a program or a file is missing). See {log}." | tee -a {log} >&2
                exit 1
            fi

            if [ "$status" -eq 0 ]; then
                outcome="GeneMark-ES ran on the test slice."
            else
                outcome="GeneMark-ES stopped on the test slice (exit status $status) but not because of its licence or installation; the real run will show whether the genome suits it."
            fi
            echo "$outcome" | tee -a {log}
            {{ date -u +"%Y-%m-%dT%H:%M:%SZ"; echo "$outcome"; }} > {output.check}
            """
