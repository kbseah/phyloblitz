import logging
import requests
import os.path
from packaging.version import Version

logger = logging.getLogger(__name__)


def list_versions():
    """List available versions of reference databases on Zenodo.

    :return: A dict of versions with their metadata keyed by version number,
        and the latest version number
    :rtype: tuple
    """
    url = "https://zenodo.org/api/records/18627380/versions"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        versions = {}
        for hit in data["hits"]["hits"]:
            versions[hit["metadata"]["version"]] = {
                "version": hit["metadata"]["version"],
                "doi": hit["metadata"]["doi"],
                "created": hit["created"],
                "files": [],
            }
            for file in hit["files"]:
                if "SSU" in file["key"]:
                    marker = "SSU"
                elif "LSU" in file["key"]:
                    marker = "LSU"
                else:
                    logger.warning(f"Unknown marker type in file: {file['key']}")
                    marker = "Unknown"
                versions[hit["metadata"]["version"]]["files"].append(
                    {
                        "marker": marker,
                        "filename": file["key"],
                        "size": file["size"],
                        "checksum": file["checksum"],
                        "download_url": file["links"]["self"],
                    }
                )
        return versions, max([Version(v) for v in versions])
    else:
        logger.error(f"Failed to fetch versions: {response.status_code}")
        return [], None


def get_file(versions, which_db, db_version, outdir=".", dryrun=False, overwrite=False):
    """Download reference database file from Zenodo

    :param versions: A dict of versions with their metadata keyed by version number, produced by list_versions()
    :param which_db: Which database to download, either "SSU" or "LSU"
    :param db_version: Version number of database to download
    :param outdir: Directory to save the downloaded file, default is current directory
    :param dryrun: If True, only print the download URL without downloading the file
    :param overwrite: If True, overwrite existing file if it already exists
    :return: Path to the downloaded file, or None if dryrun or download failed
    """
    if db_version not in versions:
        raise ValueError(f"Version {db_version!s} not found in available versions.")
    if which_db not in [file["marker"] for file in versions[db_version]["files"]]:
        raise ValueError(f"{which_db!s} database not found in version {db_version!s}.")
    for file in versions[db_version]["files"]:
        if file["marker"] == which_db:
            logger.info(f"Downloading {which_db} database version {db_version}...")
            outpath = os.path.join(outdir, file["filename"])
            if os.path.exists(outpath):
                if overwrite:
                    logger.warning(
                        f"File {outpath!s} already exists, but overwrite is enabled. Overwriting."
                    )
                else:
                    raise FileExistsError(f"File {outpath!s} already exists, but overwrite disabled. Skipping download.")
            logger.debug(f"Saving to {outpath!s}")
            if dryrun:
                logger.info(f"Dry run: would download from {file['download_url']!s}")
                return None
            response = requests.get(file["download_url"], stream=True)
            if response.status_code == 200:
                with open(outpath, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                logger.info(f"Downloaded {file['filename']} successfully.")
                # TODO: Verify checksum here
                return outpath
            else:
                raise ConnectionError(f"Failed to download file: {response.status_code!s}")
    return None
