# Release v0.5.0

Package version, changelog and annotated Git tag all use 0.5.0 / v0.5.0.
The delivery includes a Git bundle with the prepared commit and annotated tag.
To preserve it, clone the bundle and set your repository URL:

```bash
git clone ../release-assets/logos-theory-b-v0.5.0.bundle logos-theory-b
cd logos-theory-b
git remote set-url origin https://github.com/YOUR_ACCOUNT/YOUR_REPOSITORY.git
git push -u origin main
git push origin v0.5.0
```

If uploading extracted files into a different Git history, create the tag after
committing them: `git tag -a v0.5.0 -m "LOGOS synthetic release v0.5.0"`.
A ZIP file upload alone does not create GitHub tags or releases.
Use github-topics.json to set the requested repository topics in GitHub Settings.
No remote repository or hosted GitHub release is assumed to exist.

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

This is a local prepared release. The bundled tag has not been published to a
remote repository, and hosted GitHub Actions execution is not claimed.
