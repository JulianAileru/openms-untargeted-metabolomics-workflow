import csv


def read_samples(samplesheet):
    with open(samplesheet, newline="") as handle:
        return [row["sample_name"].strip() for row in csv.DictReader(handle, delimiter="\t")]


SAMPLES = read_samples(config["samplesheet"])
