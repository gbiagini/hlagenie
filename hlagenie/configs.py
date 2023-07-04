import yaml


# load the configuration file
with open("hlagenie/configurations.yml", "r") as f:
    config = yaml.safe_load(f)
