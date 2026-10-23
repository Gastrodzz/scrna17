# Deploying the scRNA-seq Platform on a VPS

This guide deploys the app for a small group of testers with **per-user
isolation and hard RAM/CPU caps**, so one big dataset can't take down the box.

Architecture: **ShinyProxy** puts each logged-in user in their **own Docker
container** (capped RAM/CPU). All containers share one persistent volume for the
SQLite accounts DB + per-project session files, so data is shared and durable.

```
                 ┌────────────── VPS ──────────────┐
  users ──HTTPS──► ShinyProxy (:8080)              │
                 │     ├─ container userA (cap 8g)  │
                 │     ├─ container userB (cap 8g)  │  each isolated
                 │     └─ ...  (max-instances = 8)  │
                 │            │                      │
                 │     /srv/scrna/data  (shared volume: DB + sessions)
                 └─────────────────────────────────┘
```

## 0. Prerequisites
- A Linux VPS (Ubuntu 22.04+). See the RAM/CPU sizing note at the bottom.
- Docker + Java 17 installed (`apt install docker.io openjdk-17-jre`).

## 1. Build the app image
From the project root (where the `Dockerfile` is):
```bash
docker build -t scrna-platform:latest .
```
First build is slow (installs the full R/Bioconductor stack via `install_dependencies.R`).
Subsequent builds are cached unless `install_dependencies.R` changes.

> **Match the base image to your R version.** The `Dockerfile` pins
> `bioconductor/bioconductor_docker:RELEASE_3_21` (R 4.5). If you use a different
> R/Bioconductor version, change that tag.

## 2. Create the shared, persistent data volume
```bash
sudo mkdir -p /srv/scrna/data
sudo chown -R 1000:1000 /srv/scrna/data   # container's default user
```
This holds `scrna4_auth.db` (accounts/projects), `user_sessions/` (per-project
figures/tables/analysis RDS) and `.keyfile` (API-key encryption master key).
**Back this directory up** — it is all the user data.

## 3. Create the Docker network ShinyProxy uses
```bash
docker network create sp-net
```

## 4. Run ShinyProxy
```bash
docker run -d --name shinyproxy --restart unless-stopped \
  --net sp-net -p 8080:8080 \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v $(pwd)/deploy/application.yml:/opt/shinyproxy/application.yml \
  openanalytics/shinyproxy:3.1.1
```
Open `http://<vps-ip>:8080`. Put **HTTPS in front** for real use (Caddy or Nginx
reverse-proxy with Let's Encrypt).

## 5. Tune resources for your box (edit `deploy/application.yml`)
- `container-memory-limit` — per-user RAM cap. Rule of thumb:
  `(VPS RAM − 4 GB overhead) / max concurrent users`.
- `max-instances` — max simultaneous containers = your queue. Extra users wait.
- Sizing for the target "10 users × up to 100k cells":

  | VPS RAM | mem-limit/user | max-instances | Notes |
  |--------:|:--------------:|:-------------:|-------|
  | 128 GB  | 12g            | 10            | comfortable, no queueing |
  | 64 GB   | 7g             | 8             | good value; add 32–64 GB swap |
  | 48–64 GB (with BPCells) | 4g | 10        | after the on-disk-matrix refactor |

  Always add generous **swap** (`fallocate -l 48G /swapfile …`) as an OOM safety net.

## Updating the app
```bash
git pull            # or copy new files
docker build -t scrna-platform:latest .
# ShinyProxy launches the new image for the next session automatically.
```
User data is untouched (it lives in the mounted volume, not the image).

## Notes
- The app's own login gates access, so ShinyProxy auth is set to `none`. For an
  extra outer gate during private testing, switch it to `simple` in `application.yml`.
- `container-env: R_MAX_VSIZE` gives heavy R jobs headroom within the cap.
- Performance parity with a Mac comes from the Bioconductor base image's
  **OpenBLAS** + a matching R version + enough cores/RAM.
