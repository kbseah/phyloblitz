import logging
import requests

logger = logging.getLogger(__name__)


def list_versions():
    """List available versions of reference databases on Zenodo.

    :return: A dict of versions with their metadata, keyed by version number
    :rtype: dict of dicts
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
        return versions
    else:
        logger.error(f"Failed to fetch versions: {response.status_code}")
        return []
