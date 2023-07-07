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
    "DQA1": "DQA1*01:01",
    "DQB1": "DQB1*05:01",
    "DPA1": "DPA1*01:03",
    "DPB1": "DPB1*01:01",
}
config["first_ten"] = {
    "A*01:01": "GSHSMRYFFT",
    "B*07:02": "GSHSMRYFYT",
    "C*01:02": "CSHSMKYFFT",
    "DRB1*01:01": "GDTRPRFLWQ",
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
config["gapped_tables"] = [f"{locus}_gapped" for locus in config["loci"]]
config["gapped_mature_tables"] = [f"{locus}_gapped_mature" for locus in config["loci"]]
config["ungapped_tables"] = [f"{locus}_ungapped" for locus in config["loci"]]
config["ungapped_mature_tables"] = [
    f"{locus}_ungapped_mature" for locus in config["loci"]
]
config["completed_tables"] = [f"{locus}_completed" for locus in config["loci"]]
config["incomplete_tables"] = [f"{locus}_incomplete" for locus in config["loci"]]
config["extended_tables"] = [f"{locus}_extended" for locus in config["loci"]]
config["position_tables"] = [f"{locus}_position" for locus in config["loci"]]
