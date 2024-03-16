#Example script for subsetting Silva derived Eukaryote database to only include Protists

bosque_euk <-metabarlist_generator(reads = reads, motus = motus, pcrs = pcrs, samples = samples2)

#write.csv(samples,"samples_bosque_bact_clean.csv")

# [1] Create Protist metabarlist
# remove kingdom_silva Chloroplastida, Fungi and Metazoa
motus$target_taxon <- !grepl("Metazoa", bosque_euk$motus$kingdom_silva)

bosque_euk_meta = subset_metabarlist(bosque_euk, table = "motus",
                                     indices = motus$target_taxon)
motus_meta<-extract_table(bosque_euk_meta, table = "motus")

motus_meta$target_taxon <- !grepl("Fungi", bosque_euk_meta$motus$kingdom_silva)

bosque_euk_fungi = subset_metabarlist(bosque_euk_meta, table = "motus",
                                      indices = motus_meta$target_taxon)

motus_fungi<-extract_table(bosque_euk_fungi, table = "motus")

motus_fungi$target_taxon <- !grepl("Chloroplastida", bosque_euk_fungi$motus$kingdom_silva)

bosque_euk_protist = subset_metabarlist(bosque_euk_fungi, table = "motus",
                                        indices = motus_fungi$target_taxon)