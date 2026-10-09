#!/bin/sh
# -----------------------------------------------------------------------------
# push_github_mirror.sh - publish the GitHub mirror at
#     github.com/xxcai1/ncplugin-SasCSNS
# from this (ihep-hosted) repository. The mirror is ONE-WAY: ihep -> github.
#
# doc/ holds draft-paper material and must not become public, so the mirror
# strips doc/ from ALL history -- except doc/wiki/, which is published.
# Consequences, by design:
#   * commit SHAs on github DIFFER from the ihep SHAs (rewritten history),
#   * the push is a force push (the mirror is defined by this script),
#   * the mirror is deterministic: re-running the script on the same ihep
#     history yields the same github SHAs (fixed date on the link-fix commit).
# The pre-mirror state of the github repo is preserved on the branch
# archive/pre-mirror-2026-10.
# -----------------------------------------------------------------------------
set -e

SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK="${TMPDIR:-/tmp}/ncplugin_github_mirror"
GH=git@github.com:xxcai1/ncplugin-SasCSNS.git
SPEC_URL=https://code.ihep.ac.cn/cinema-developers/ncplugin-sascsns/-/blob/main/doc/anisotropic_directload2d.pdf

rm -rf "$WORK"
git clone --no-local --quiet "$SRC" "$WORK"
cd "$WORK"

# strip doc/ except doc/wiki/, across all history
git filter-repo --force --invert-paths --path-regex '^doc/(?!wiki/)'

# point the README's spec link at the canonical ihep location (the file is
# not in the mirror); fixed dates keep the commit SHA reproducible
sed -i "s#](doc/anisotropic_directload2d.pdf)#]($SPEC_URL)#" README.md
if ! git diff --quiet; then
    GIT_AUTHOR_DATE='2026-10-09T00:00:00 +0000' \
    GIT_COMMITTER_DATE='2026-10-09T00:00:00 +0000' \
    git commit -aqm 'mirror: link the specification at its canonical ihep location (doc/ is stripped from this mirror)'
fi

git remote add github "$GH"
git push -q --force github main
echo "mirror pushed: $(git rev-parse --short main) (github main)"
echo "stripped check: $(git ls-tree -r --name-only main | grep -c '^doc/') doc/ files remain, of which $(git ls-tree -r --name-only main | grep -c '^doc/wiki/') are doc/wiki/"
