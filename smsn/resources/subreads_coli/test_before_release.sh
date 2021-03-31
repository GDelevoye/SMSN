#informal test that must be ran manually
#conda activate smsn # --> We must be in the conda env
echo "Informal test that must be ran manually before releasing"
echo "Ensure that you're in the smsn conda environment"
echo "cmd1 = smsn --bam ./HTVEG_Coli.subreads.bam --reference ./../e_coli_O157.fasta --output_csv ./methylationout_HTVEG_Coli.csv --verbosity INFO --add_context True --progress_bar True --nb_proc 12 --min_subreads 50 --idQvs True --CCS ./CCS_HTVEG.bam --min_identity 0.99"
smsn --bam ./HTVEG_Coli.subreads.bam --reference ./../e_coli_O157.fasta --output_csv ./methylationout_HTVEG_Coli.csv --verbosity INFO --add_context True --progress_bar True --nb_proc 12 --min_subreads 50 --idQvs True --CCS ./CCS_HTVEG.bam --min_identity 0.99
echo "cmd2 = smsn --bam ./HT2_Coli.subreads.bam --reference ./../e_coli_O157.fasta --output_csv ./methylationout_HT2_Coli.csv --verbosity DEBUG --add_context False --progress_bar False --nb_proc 1 --min_subreads 20 --idQvs True --model SP2-C2"
smsn --bam ./HT2_Coli.subreads.bam --reference ./../e_coli_O157.fasta --output_csv ./methylationout_HT2_Coli.csv --verbosity DEBUG --add_context False --progress_bar False --nb_proc 1 --min_subreads 20 --idQvs True --model SP2-C2
echo "cmd3 = smsn --bam ./HT6_Coli.subreads.bam --reference ./../e_coli_O157.fasta --output_csv ./methylationout_HT6_Coli.csv --verbosity WARNING --add_context True --progress_bar True --nb_proc 15 --min_subreads 20 --idQvs False --model SP2-C2"
smsn --bam ./HT6_Coli.subreads.bam --reference ./../e_coli_O157.fasta --output_csv ./methylationout_HT6_Coli.csv --verbosity WARNING --add_context True --progress_bar True --nb_proc 15 --min_subreads 15 --idQvs False --model SP2-C2
