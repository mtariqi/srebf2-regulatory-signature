# Defined interactively
function count_pilot_strand_check
    for strand_mode in 0 1 2
        echo "Starting strand mode $strand_mode"

        featureCounts \
                        -T 4 \
                        -p --countReadPairs \
                        -B -C \
                        -s $strand_mode \
                        -t exon \
                        -g gene_id \
                        -a data/raw/reference/gencode_v47/gencode.v47.primary_assembly.annotation.gtf \
                        -o results/qc/pilot/strand_check/SRR28966297_s$strand_mode.txt \
                        data/processed/alignment/SRR28966297.sorted.bam \
                        > logs/SRR28966297_featureCounts_s$strand_mode.log 2>&1

        set count_status $status
        echo "Strand mode $strand_mode exit status: $count_status"

        if test $count_status -ne 0
            tail -n 30 logs/SRR28966297_featureCounts_s$strand_mode.log
            return 1
        end
    end
end
