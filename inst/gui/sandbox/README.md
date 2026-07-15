# Spectreasy onboarding sandbox

This local-only controller tests onboarding without reading or modifying the host R installation.

- The real sandbox stores Micromamba, R, installed packages, and its empty project under `~/.cache/spectreasy-onboarding-sandbox/` by default. Set `SPECTREASY_SANDBOX_ROOT` to override it.
- Reset removes only the sandbox R environment, its R library, and sandbox project. It retains Micromamba's download cache so later R-runtime tests are faster.
- Delete all sandbox files removes the runtime, libraries, project, tools, and download cache. A clean R plus the compiled Spectreasy dependency stack can use roughly 2 GB during testing.
- Preview states do not write anything. They exist for instant visual testing of missing, outdated, broken, and ready states.
- The controller binds only to `127.0.0.1` and exposes fixed actions; it does not accept arbitrary shell commands.

Run the controller and frontend in separate terminals:

```sh
cd inst/gui
npm run sandbox
npm run dev -- --host 127.0.0.1
```

Then open:

```text
http://127.0.0.1:5174/?sandbox=1&api=http://127.0.0.1:8787
```
