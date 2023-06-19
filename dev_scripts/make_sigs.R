library(xCell2)

# RNA-seq references
# TODO: Add Mahmoud's reference
print("Kassandra Blood Reference")
kass_blood_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/kass_blood_ref.rds")
kass_blood_sigs <- xCell2Train(kass_blood_ref$ref, kass_blood_ref$labels, data_type = "rnaseq", lineage_file = kass_blood_ref$lineage_file)
saveRDS(kass_blood_sigs, "/bigdata/almogangel/xCell2_data/dev_data/kass_blood_sigs.rds")
print("Done")

print("Kassandra Tumor Reference")
kass_tumor_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_ref.rds")
kass_tumor_sigs <- xCell2Train(kass_tumor_ref$ref, kass_tumor_ref$labels, data_type = "rnaseq", lineage_file = kass_tumor_ref$lineage_file)
saveRDS(kass_tumor_sigs, "/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_sigs.rds")
print("Done")

print("BlueprintEncode Reference")
bp_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/bp_ref.rds")
bp_sigs <- xCell2Train(bp_ref$ref, bp_ref$labels, data_type = "rnaseq", lineage_file = bp_ref$lineage_file)
saveRDS(bp_sigs, "/bigdata/almogangel/xCell2_data/dev_data/bp_sigs.rds")
print("Done")



# Array references
print("LM22 Reference")
lm22_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/lm22_ref.rds")
lm22_sigs <- xCell2Train(lm22_ref$ref, lm22_ref$labels, data_type = "array", lineage_file = lm22_ref$lineage_file)
saveRDS(lm22_sigs, "/bigdata/almogangel/xCell2_data/dev_data/lm22_sigs.rds")
print("Done")



# scRNA-seq references
# TODO: Add Tabula Sapiens reference
