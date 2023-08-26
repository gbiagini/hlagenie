#!/usr/bin/env python3

config = {}

config["imgt_version"] = "Latest"
config["DEFAULT_CACHE_SIZE"] = 1000
config["loci"] = [
    "A",
    "B",
    "C",
    "DRB1",
    "DQA1",
    "DQB1",
    "DPA1",
    "DPB1",
]
config["refseq"] = {
    "A": "A*01:01",
    "B": "B*07:02",
    "C": "C*01:02",
    "DRB1": "DRB1*01:01",
    "DRB3": "DRB3*01:01",
    "DRB4": "DRB4*01:01",
    "DRB5": "DRB5*01:01",
    "DQA1": "DQA1*01:01",
    "DQB1": "DQB1*05:01",
    "DPA1": "DPA1*01:03",
    "DPB1": "DPB1*01:01",
}
config["refseq_full"] = {
    "A": "A*01:01:01:01",
    "B": "B*07:02:01:01",
    "C": "C*01:02:01:01",
    "DRB1": "DRB1*01:01:01:01",
    "DRB3": "DRB3*01:01:02:01",
    "DRB4": "DRB4*01:01:01:01",
    "DRB5": "DRB5*01:01:01:01",
    "DQA1": "DQA1*01:01:01:01",
    "DQB1": "DQB1*05:01:01:01",
    "DPA1": "DPA1*01:03:01:01",
    "DPB1": "DPB1*01:01:01:01",
}
config["first_ten"] = {
    "A*01:01": "GSHSMRYFFT",
    "B*07:02": "GSHSMRYFYT",
    "C*01:02": "CSHSMKYFFT",
    "DRB1*01:01": "GDTRPRFLWQ",
    "DRB3*01:01": "GDTRPRFLEL",
    "DRB4*01:01": "GDTQPRFLEQ",
    "DRB5*01:01": "GDTRPRFLQQ",
    "DQA1*01:01": "EDIVADHVAS",
    "DQB1*05:01": "RDSPEDFVYQ",
    "DPA1*01:03": "IKADHVSTYA",
    "DPB1*01:01": "RATPENYVYQ",
}
config["ard_last_ten"] = {
    "A*01:01": "NGKETLQRTD",
    "B*07:02": "NGKDKLERAD",
    "C*01:02": "NGKETLQRAE",
    "DRB1*01:01": "GESFTVQRRV",
    "DRB3*01:01": "GESFTVQRRV",
    "DRB4*01:01": "VESFTVQRRV",
    "DRB5*01:01": "GESFTVQRRV",
    "DQA1*01:01": "RYNSTAATNE",
    "DQB1*05:01": "AYRGILQRRV",
    "DPA1*01:03": "RSNHTQATND",
    "DPB1*01:01": "DEAVTLQRRV",
}
config["xrd_last_ten"] = {
    "A*01:01": "LPKPLTLRWE",
    "B*07:02": "LPKPLTLRWE",
    "C*01:02": "LPEPLTLRWE",
    "DRB1*01:01": "VTSPLTVEWR",
    "DRB3*01:01": "VTSALTVEWR",
    "DRB4*01:01": "MMSPLTVQWS",
    "DRB5*01:01": "VTSPLTVEWR",
    "DQA1*01:01": "LDQPLLKHWE",
    "DQB1*05:01": "LQSPITVEWR",
    "DPA1*01:03": "LDQPLLKHWE",
    "DPB1*01:01": "LDSPVTVEWK",
}
config["expression_chars"] = [
    "N",
    "Q",
    "L",
    "S",
]
# TODO: rename gapped_tables to gapped_tables_prot
config["gapped_tables"] = [f"{locus}_gapped" for locus in config["loci"]]
config["gapped_mature_tables"] = [f"{locus}_gapped_mature" for locus in config["loci"]]
# TODO: rename ungapped_tables to ungapped_tables_prot
config["ungapped_tables"] = [f"{locus}_ungapped" for locus in config["loci"]]
config["ungapped_mature_tables"] = [
    f"{locus}_ungapped_mature" for locus in config["loci"]
]
config["ungapped_nuc_tables"] = [f"{locus}_ungapped_nuc" for locus in config["loci"]]
config["gapped_nuc_tables"] = [f"{locus}_gapped_nuc" for locus in config["loci"]]
config["position_tables"] = [f"{locus}_position" for locus in config["loci"]]
