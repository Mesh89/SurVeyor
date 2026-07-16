import hashlib
import json
import os
import re


MANIFEST_FILENAME = "training_manifest.json"
VCF_HEADER_KEY = "SurVeyorTrainingSetSHA256"
_READ_BUFFER_SIZE = 1024 * 1024
_SHA256_PATTERN = re.compile(r"[0-9a-fA-F]{64}")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as input_file:
        for chunk in iter(lambda: input_file.read(_READ_BUFFER_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_training_manifest(training_prefixes, model_name, cross_species):
    file_specs = [
        ("variants", ".vcf.gz"),
        ("statistics", ".stats"),
        ("genotypes", ".gts"),
    ]
    if model_name in ("ALL", "INS_TO_DUP", "INS_TO_DUP_LARGE"):
        file_specs.insert(1, ("insertions_to_duplications", ".INS_TO_DUP.vcf.gz"))

    samples = []
    for sample_index, prefix in enumerate(training_prefixes):
        sample_files = []
        for role, suffix in file_specs:
            path = prefix + suffix
            sample_files.append({
                "role": role,
                "path": path,
                "size_bytes": os.path.getsize(path),
                "sha256": sha256_file(path),
            })
        samples.append({
            "sample_index": sample_index,
            "prefix": prefix,
            "files": sample_files,
        })

    # Keep this aggregate independent of where the input files live. Paths remain
    # in the manifest for auditing, while the training-set ID reflects content,
    # file roles, sample grouping, and sample order.
    training_set_contents = [
        [
            {"role": sample_file["role"], "sha256": sample_file["sha256"]}
            for sample_file in sample["files"]
        ]
        for sample in samples
    ]
    training_set_payload = json.dumps(
        training_set_contents,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")

    return {
        "manifest_version": 1,
        "hash_algorithm": "sha256",
        "training_set_sha256": hashlib.sha256(training_set_payload).hexdigest(),
        "model_name": model_name,
        "cross_species": cross_species,
        "samples": samples,
    }


def write_training_manifest(output_dir, manifest):
    os.makedirs(output_dir, exist_ok=True)
    manifest_path = os.path.join(output_dir, MANIFEST_FILENAME)
    temporary_path = manifest_path + ".tmp"
    with open(temporary_path, "w", encoding="utf-8") as output_file:
        json.dump(manifest, output_file, indent=2, ensure_ascii=False)
        output_file.write("\n")
    os.replace(temporary_path, manifest_path)
    return manifest_path


def load_training_set_sha256(model_dir):
    manifest_path = os.path.join(model_dir, MANIFEST_FILENAME)
    if not os.path.exists(manifest_path):
        return None

    try:
        with open(manifest_path, encoding="utf-8") as input_file:
            manifest = json.load(input_file)
    except (OSError, json.JSONDecodeError) as error:
        raise RuntimeError(
            f"Could not read training manifest: {manifest_path}"
        ) from error

    training_set_sha256 = manifest.get("training_set_sha256")
    if (
        manifest.get("hash_algorithm") != "sha256"
        or not isinstance(training_set_sha256, str)
        or _SHA256_PATTERN.fullmatch(training_set_sha256) is None
    ):
        raise RuntimeError(
            f"Training manifest has an invalid training_set_sha256: {manifest_path}"
        )
    return training_set_sha256.lower()


def training_set_vcf_header_line(training_set_sha256):
    if _SHA256_PATTERN.fullmatch(training_set_sha256) is None:
        raise ValueError("training_set_sha256 must be a 64-character hexadecimal digest")
    return f"##{VCF_HEADER_KEY}={training_set_sha256.lower()}"


def set_training_set_vcf_header(vcf_header, training_set_sha256):
    for header_record in list(vcf_header.records):
        if header_record.key == VCF_HEADER_KEY:
            header_record.remove()
    vcf_header.add_line(training_set_vcf_header_line(training_set_sha256))
