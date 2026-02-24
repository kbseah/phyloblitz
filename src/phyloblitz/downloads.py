import hashlib
import logging
import os.path

import requests
from packaging.version import Version
from tqdm import tqdm

logger = logging.getLogger(__name__)

DB_URL_SILVA = "https://zenodo.org/api/records/18627380/versions"


def list_versions():
    """List available versions of reference databases on Zenodo.

    :return: A dict of versions with their metadata keyed by version number,
        and the latest version number
    :rtype: tuple
    """
    response = requests.get(DB_URL_SILVA)
    if response.status_code == 200:
        data = response.json()
        versions = {}
        for hit in data["hits"]["hits"]:
            versions[hit["metadata"]["version"]] = {
                "version": hit["metadata"]["version"],
                "doi": hit["metadata"]["doi"],
                "created": hit["created"],
                "files": {},
            }
            for file in hit["files"]:
                if "SSU" in file["key"]:
                    marker = "SSU"
                elif "LSU" in file["key"]:
                    marker = "LSU"
                else:
                    logger.warning("Unknown marker type in file: %s", file["key"])
                    marker = "Unknown"
                versions[hit["metadata"]["version"]]["files"][marker] = {
                    "marker": marker,
                    "filename": file["key"],
                    "size": file["size"],
                    "checksum": file["checksum"],
                    "download_url": file["links"]["self"],
                }
        return versions, max([Version(v) for v in versions])
    logger.error("Failed to fetch versions: %s", str(response.status_code))
    return [], None


def get_file(versions, which_db, db_version, outdir=".", dryrun=False, overwrite=False):
    """Download reference database file from Zenodo.

    :param versions: A dict of versions with their metadata keyed by version
        number, produced by list_versions()
    :param which_db: Which database to download, either "SSU" or "LSU"
    :param db_version: Version number of database to download
    :param outdir: Directory to save the downloaded file, default is current directory
    :param dryrun: If True, only print the download URL without downloading the file
    :param overwrite: If True, overwrite existing file if it already exists
    :return: Path to the downloaded file, or the expected path if dryrun is True
    """
    if db_version not in versions:
        msg = f"Version {db_version!s} not found in available versions."
        raise ValueError(msg)
    if which_db not in versions[db_version]["files"]:
        msg = f"{which_db!s} database not found in version {db_version!s}."
        raise ValueError(msg)
    logger.info("Downloading %s database version %s ...", which_db, db_version)
    outpath = os.path.join(outdir, versions[db_version]["files"][which_db]["filename"])
    if os.path.exists(outpath):
        if overwrite:
            logger.warning(
                "File %s already exists, but overwrite is enabled. Overwriting.",
                str(outpath),
            )
        else:
            msg = f"File {outpath!s} already exists, but overwrite disabled. Skipping download."
            raise FileExistsError(msg)
    logger.debug("Saving to %s", str(outpath))
    if dryrun:
        logger.info(
            "Dry run: would download from %s",
            versions[db_version]["files"][which_db]["download_url"],
        )
        return outpath
    response = requests.get(
        versions[db_version]["files"][which_db]["download_url"], stream=True,
    )
    if response.status_code == 200:
        total_size = int(response.headers.get("Content-Length", 0))
        logger.debug("Expected total file size: %d bytes", total_size)
        with (
            open(outpath, "wb") as f,
            tqdm(
                total=total_size,
                unit="B",
                bar_format="{l_bar}{bar:20}{r_bar}",
                colour="green",
                unit_scale=True,
                desc=versions[db_version]["files"][which_db]["filename"],
            ) as pbar,
        ):
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                pbar.update(len(chunk))
        logger.info(
            "Downloaded %s successfully.",
            versions[db_version]['files'][which_db]["filename"],
        )
        return outpath
    raise ConnectionError(f"Failed to download file: {response.status_code!s}")


def check_md5sum_file(versions, which_db, db_version, outpath):
    """Check MD5 checksum of downloaded file.

    :param versions: A dict of versions with their metadata keyed by version
        number, produced by list_versions()
    :param which_db: Which database downloaded
    :param db_version: Version number of database downloaded
    :param outpath: Path where file was downloaded
    :return: True if checksum matches value in metadata
    :rtype: bool
    """
    checksum_type, expected_checksum = versions[db_version]["files"][which_db][
        "checksum"
    ].split(":")
    if checksum_type != "md5":
        msg = f"Checksum type {checksum_type!s} not supported. Only md5 is supported."
        raise NotImplementedError(msg)
    md5_hash = hashlib.md5()
    with open(outpath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)
    calculated_md5 = md5_hash.hexdigest()
    logger.debug("Calculated MD5 checksum for %s: %s", outpath, calculated_md5)
    logger.debug("Expected MD5 checksum from metadata: %s", expected_checksum)
    return calculated_md5 == expected_checksum
