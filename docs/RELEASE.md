# Release v0.5.1

Package version, changelog and annotated Git tag all use 0.5.1 / v0.5.1.
The delivery includes a Git bundle with the prepared commit and annotated tag.
To preserve it, clone the bundle and set your repository URL:

```bash
git clone ../release-assets/logos-theory-b-v0.5.1.bundle logos-theory-b
cd logos-theory-b
git remote set-url origin https://github.com/vkaralis/LOGOS-theory-B-apps.git
git push -u origin main
git push origin v0.5.1
```

If uploading extracted files into a different Git history, create the tag after
committing them: `git tag -a v0.5.1 -m "LOGOS synthetic release v0.5.1"`.
A ZIP file upload alone does not create GitHub tags or releases.
Use github-topics.json to set the requested repository topics in GitHub Settings.
The public repository is [vkaralis/LOGOS-theory-B-apps](https://github.com/vkaralis/LOGOS-theory-B-apps).
Published releases are available on its [Releases page](https://github.com/vkaralis/LOGOS-theory-B-apps/releases).

## GitHub About and topics

Set the repository description to the exact text in `GITHUB_ABOUT.txt` and add
all eleven topics from `github-topics.json`. These local files do not change remote GitHub settings.
The longer research description appears at the beginning of `README.md`,
followed by the explicit scope of the synthetic release.

The ZIP contains the repository folder and its matching annotated-tag Git
bundle. Upload the contents of the repository folder as the repository files;
do not include the external `release-assets/` directory in the source tree.
The manuscript and original source-data archive are excluded from this release.
The MIT license applies to the bundled code, documentation and generated
synthetic datasets; retain its copyright and permission notice on reuse.

Release v0.5.0 and its annotated tag were published and passed hosted GitHub
Actions. Version v0.5.1 updates documentation and adds release-inventory
verification. Its source files, ZIP and Git bundle are prepared together and
its hosted checks are verified after publication. Earlier published tags and
release artifacts remain historical snapshots.
