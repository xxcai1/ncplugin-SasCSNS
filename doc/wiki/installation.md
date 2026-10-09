# Installation

How to install the **SasCSNS** NCrystal plugin and verify that `python` can see it. For what the plugin does, see [Home](home); for the NCMAT data it consumes, see [Data format](data-format).

## Prerequisites

| Requirement | Version | Where it is declared |
|---|---|---|
| NCrystal (`ncrystal-core`) | `>=4.2.0,<4.3` | runtime dependency and build requirement in [pyproject.toml](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/pyproject.toml) |
| `ncrystal-pypluginmgr` | `>=0.0.5` | runtime dependency in [pyproject.toml](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/pyproject.toml) |
| C++ compiler | any supported by your CMake/NCrystal setup | the CMake project is declared `LANGUAGES CXX` in [CMakeLists.txt](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/CMakeLists.txt) |
| CMake | `>=3.20` | `cmake_minimum_required(VERSION 3.20...3.31)` in [CMakeLists.txt](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/CMakeLists.txt) |
| Python | 3 | only the `Programming Language :: Python :: 3` classifier is declared for Python (no `requires-python` bound); no minimum version is pinned in `pyproject.toml` |

`scikit-build-core>=0.10` is the build backend; `pip` fetches it automatically during the build.

NCrystal itself can be provided either by conda or by the `ncrystal-core` Python wheel: the CMake code locates the NCrystal package via `ncrystal-config --show cmakedir` if `NCrystal_DIR` is not set, which the comments in `CMakeLists.txt` note is what makes the build "work with ncrystal-core installed via python wheels". Only the CMake-level `find_package( NCrystal 4.0.0 REQUIRED )` is looser than the pip pin above; the `pip install .` route always enforces `>=4.2.0,<4.3`.

## Install

Run all commands with the NCrystal environment active — the same environment that `python` and `pip` resolve to.

1. Activate the environment that provides NCrystal, e.g. for conda:

   ```sh
   conda activate <env-with-ncrystal>
   ```

2. Clone the repository:

   ```sh
   git clone https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns.git
   cd ncplugin-sascsns
   ```

3. Build and install (this invokes CMake automatically through scikit-build-core, and installs the plugin package `ncrystal-plugin-SasCSNS` into the active environment):

   ```sh
   pip install .
   ```

That is the whole install: the README instructs to "Run from the repo root with the NCrystal+plugin environment active (`python` must see the plugin, i.e. `pip install .` first)". After the install, the environment's site-packages contain:

```text
ncrystal_plugin_SasCSNS/
├── __init__.py
├── plugins/libNCPlugin_SasCSNS.so   # the compiled plugin
└── data/sascsns_sio2_spheres.ncmat  # bundled example data file
```

The plugin name is the single line of [ncplugin_name.txt](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/ncplugin_name.txt) (`SasCSNS`); CMake aborts the build if it disagrees with the `ncrystal-plugin-SasCSNS` project name in `pyproject.toml`.

### Plain CMake build (without pip)

`pip install .` is the supported route, because it is the one that installs the plugin where NCrystal's plugin manager looks for it. A bare CMake build also works for compiling and inspecting:

```sh
cmake -S . -B build-standalone -DCMAKE_INSTALL_PREFIX=<prefix>
cmake --build build-standalone
cmake --install build-standalone
```

This uses the fallback install layout in `CMakeLists.txt` (`lib/` and `data/` under the prefix). <!-- TODO-WIKI: whether a non-pip install (lib/ under a custom prefix) is auto-discovered by NCrystal is not documented in the repo; confirm before recommending it to users. -->

## Verify the installation

### Plugin self-test

With the install environment active:

```sh
ncrystal-pluginmanager --test SasCSNS
```

A successful run ends like this (captured with NCrystal 4.2.12 and plugin 1.0.0; the numeric values are regression anchors and can change between plugin/NCrystal versions):

```text
NCrystal: Required plugin "SasCSNS" was indeed available (dynamic).
NCrystal: Launching plugin test function "test_SasCSNS":
NCrystal: plugin::SasCSNS: Testing plugin
...
NCrystal: plugin::SasCSNS: All tests of plugin were successful!
NCrystal: End of plugin test function "test_SasCSNS".
Testing load of "plugins::SasCSNS/sascsns_sio2_spheres.ncmat"
  -> createTextData
  -> createInfo
  -> createScatter
  -> createAbsorption
All ok
```

The `--test` flag is the documented way to "load and test a given dynamic plugin" (`ncrystal-pluginmanager --help`). Note that `ncrystal-pluginmanager --list` is **not** a supported flag in `ncrystal-pypluginmgr` 0.0.5 — the CLI treats `--list` as a plugin name and fails with `NCLogicError`.

### Python API smoke test

```python
import NCrystal
mat = NCrystal.createLoadedMaterial('plugins::SasCSNS/sascsns_sio2_spheres.ncmat')
print('SasCSNS OK:', mat.scatter is not None)
```

The `plugins::SasCSNS/...` path resolves to the data file bundled with the installed plugin package (the same load path used by the plugin's own test above). If this prints `SasCSNS OK: True`, NCrystal discovered the plugin and can build a scatter model from its `@CUSTOM_SASCSNS` section. The equivalent manual check on any `@CUSTOM_SASCSNS` file is `python -c "import NCrystal; NCrystal.createLoadedMaterial('<your-file>.ncmat')"`.

Running the full validation ladder afterwards is described on the [Validation](validation) page; hands-on usage starts at [Tutorial 1D](tutorial-1d) or [Tutorial 2D](tutorial-2d).

## Notes

- **`python` must see the plugin.** All NCrystal work goes through the Python API: the example, figure and plugin-driving validation scripts under [script/](https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/tree/main/script) all import NCrystal (the only exceptions are `sasview2ncmat.py`, a pure-text converter, and `validate_directload2d.py`, which checks the Python model against closed forms without NCrystal), so the interpreter, NCrystal and the installed plugin must live in the same environment. If you install with `pip` from one environment and run `python` from another, the plugin will not be found.
- **Plugin discovery.** The `ncrystal-pluginmanager` command "is invoked by the NCrystal library in order to discover any plugins available in the environment" (`ncrystal-pluginmanager --help`); the plugin is shipped as the pip package `ncrystal-plugin-SasCSNS`, which is what this discovery mechanism finds. (pip normalises distribution names, so the underscore spelling `ncrystal_plugin_SasCSNS` — the installed module directory — is equivalent on the command line.) No Python entry points are used.
- **Expected warning.** Loading any file with a `@CUSTOM_SASCSNS` section prints `NCrystal WARNING: Loading NCMAT data which has @CUSTOM_ section(s). This is OK if intended.` This is normal and appears in the plugin's own test output.
- **Disabling SANS.** The plugin respects the NCrystal `sans` request parameter: with `sans=0` the plugin factory declines the request (so the SANS model is not added); multi-phase materials (`@OTHERPHASES`) are likewise declined and left to NCrystal's builtin factories.
- **Conda notes.** A conda environment is a convenient way to provide NCrystal (`ncrystal`, `ncrystal-core`, `ncrystal-lib`, `ncrystal-python` packages); the plugin itself is then installed with `pip install .` on top. The `build-conda/` directory sometimes seen in checkouts is *not* an environment setup — it is a plain, gitignored CMake build tree from a local build.
- **Reference environment.** One development setup observed to pass all checks: conda-provided NCrystal 4.2.12 + pip-installed `ncrystal-pypluginmgr` 0.0.5 + plugin 1.0.0 on Python 3.10. These are environment observations, not requirements; only the versions in the table above are declared by the project.

## Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| `ncrystal-pluginmanager --test SasCSNS` says the plugin was not available, or `@CUSTOM_SASCSNS` files fail to load | Plugin not installed into the environment you are running — the classic case is `pip install .` done with a different `pip`/`python` than the active one. Re-run `pip install .` with the NCrystal environment active, then re-test. |
| Load of a `@CUSTOM_SASCSNS` file produces no SANS model although the plugin tests fine | The request may have SANS disabled (`sans=0`) or the material is multi-phase (`@OTHERPHASES`), both of which the plugin declines. Remove `sans=0`, or use a single-phase file. |
| `pip install .` fails resolving NCrystal | The active environment's NCrystal is outside the pinned range (`ncrystal-core>=4.2.0,<4.3`) or missing. Install a matching NCrystal first (conda or `pip install ncrystal-core`). |
| Build cannot find the NCrystal CMake package | Make sure `ncrystal-config` from the active environment is on `PATH` (CMake calls `ncrystal-config --show cmakedir`), or pass `-DNCrystal_DIR=<path>` explicitly. |
| `ncrystal-pluginmanager --list` errors with `NCLogicError` | `--list` does not exist in `ncrystal-pypluginmgr` 0.0.5; use `ncrystal-pluginmanager --test SasCSNS`. |
| Converter output fails to load on an older NCrystal | `script/sasview2ncmat.py` writes NCMAT v6 (1D) and NCMAT v7 (2D DirectLoad2D) files; the reading NCrystal must support those NCMAT versions. See [Data format](data-format). |
