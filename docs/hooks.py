import glob
import os
import shutil


def copy_changelog_license_and_images(*args, **kwargs):
    shutil.copy("CHANGELOG.md", "docs/changelog.md")
    shutil.copy("LICENSE.md", "docs/license.md")
    os.makedirs("docs/images", exist_ok=True)
    for image in glob.glob("images/*.png"):
        shutil.copy(image, "docs/images/")


def exclude_agent_docs(files, *args, **kwargs):
    # docs/agents/ and docs/adr/ hold configuration for coding agents, not user documentation
    for file in list(files):
        if file.src_uri.startswith(("agents/", "adr/")):
            files.remove(file)
    return files
