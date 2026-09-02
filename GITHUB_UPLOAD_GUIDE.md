# GitHub upload guide

Repository name: `CXCL13-Tpex-HCC-niche`

1. Sign in to <https://github.com/liukairui621> and select **New repository**.
2. Enter `CXCL13-Tpex-HCC-niche`, set visibility to **Public**, and do not add a README, license, or .gitignore because they are already included.
3. Open PowerShell in the local repository folder and run:

```powershell
git init
git add .
git commit -m "Initial public release for PLOS ONE submission"
git branch -M main
git remote add origin https://github.com/liukairui621/CXCL13-Tpex-HCC-niche.git
git push -u origin main
```

4. Confirm that `https://github.com/liukairui621/CXCL13-Tpex-HCC-niche` opens without signing in.
5. Create a release such as `v1.0-submission`. For an immutable citation, archive the release with Zenodo and add the DOI to the PLOS Data Availability field if available.

Do not submit the Data Availability statement claiming repository access until step 4 succeeds.
