library(xCell2)

# Load validation datasets


# RNA-seq references
# TODO: Add Mahmoud's reference
print("Kassandra Blood Reference")
kass_blood_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/kass_blood_ref.rds")
kass_blood_sigs <- xCell2Train(kass_blood_ref$ref, kass_blood_ref$labels, data_type = "rnaseq", lineage_file = kass_blood_ref$lineage_file, filter_sigs = FALSE)
# saveRDS(kass_blood_sigs, "/bigdata/almogangel/xCell2_data/dev_data/kass_blood_sigs.rds")
saveRDS(kass_blood_sigs, "/bigdata/almogangel/xCell2_data/dev_data/kass_blood_all_sigs.rds")
print("Done")

print("Kassandra Tumor Reference")
kass_tumor_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_ref.rds")
kass_tumor_sigs <- xCell2Train(kass_tumor_ref$ref, kass_tumor_ref$labels, data_type = "rnaseq", lineage_file = kass_tumor_ref$lineage_file, filter_sigs = FALSE)
# saveRDS(kass_tumor_sigs, "/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_sigs.rds")
saveRDS(kass_tumor_sigs, "/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_all_sigs.rds")
print("Done")

print("BlueprintEncode Reference")
bp_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/bp_ref.rds")
bp_sigs <- xCell2Train(bp_ref$ref, bp_ref$labels, data_type = "rnaseq", lineage_file = bp_ref$lineage_file, filter_sigs = FALSE)
# saveRDS(bp_sigs, "/bigdata/almogangel/xCell2_data/dev_data/bp_sigs.rds")
saveRDS(bp_sigs, "/bigdata/almogangel/xCell2_data/dev_data/bp_all_sigs.rds")
print("Done")


# Array references
print("LM22 Reference")
lm22_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/lm22_ref.rds")
lm22_sigs <- xCell2Train(lm22_ref$ref, lm22_ref$labels, data_type = "array", lineage_file = lm22_ref$lineage_file, filter_sigs = FALSE)
# saveRDS(lm22_sigs, "/bigdata/almogangel/xCell2_data/dev_data/lm22_sigs.rds")
saveRDS(lm22_sigs, "/bigdata/almogangel/xCell2_data/dev_data/lm22_all_sigs.rds")
print("Done")


# scRNA-seq references
print("Tabula Sapiens Blood Reference")
ts_blood_ref <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/ts_blood_ref.rds")
ts_blood_sigs <- xCell2Train(ref = ts_blood_ref$ref, labels = ts_blood_ref$labels, data_type = "sc", lineage_file = ts_blood_ref$lineage_file, filter_sigs = FALSE)
# saveRDS(ts_blood_sigs, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/xcell2_sigs/ts_blood_sigs.rds")
saveRDS(ts_blood_sigs, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/xcell2_sigs/ts_blood_all_sigs.rds")
print("Done")

print("Pan Cancer Reference")
pan_cancer_ref <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/sc_pan_cancer_ref.rds")
pan_cancer_sigs <- xCell2Train(ref = pan_cancer_ref$ref, labels = pan_cancer_ref$labels, data_type = "sc", lineage_file = pan_cancer_ref$lineage_file, filter_sigs = FALSE)
# saveRDS(pan_cancer_sigs, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/xcell2_sigs/sc_pan_cancer_sigs.rds")
saveRDS(pan_cancer_sigs, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/xcell2_sigs/sc_pan_cancer_all_sigs.rds")
print("Done")
