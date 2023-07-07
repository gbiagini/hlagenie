import yaml
from os import path

# get the path to the configuration file
filepath = path.join(path.dirname(__file__), "configurations.yml")

# load the configuration file
with open(filepath, "r") as f:
    config = yaml.safe_load(f)
