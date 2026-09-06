"""Synthetic Makefiles build evidence: no real compilation or benchmark needed."""
import os
import pathlib
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import run_slim_regression as runner
import slim_build_provenance as provenance


class BuildProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name)
        self.runtime = self.root / "libruntime.so"
        self.runtime.write_bytes(b"same runtime")
        mock = patch.object(provenance, "output", side_effect=self.command_output)
        mock.start()
        self.addCleanup(mock.stop)

    def command_output(self, command, cwd=None):
        if command[-1] == "--version":
            return command[0] + " fixture compiler version 1"
        if command[0] == "ldd":
            return f"libruntime.so => {self.runtime} (0x1234)"
        if command[:2] == ["git", "ls-tree"]:
            return ("160000 commit abc\tthird_party/spral" if "-r" in command
                    else "040000 tree def\tthird_party")
        if command[:2] == ["git", "rev-parse"]:
            return "abc"
        if command[:2] == ["git", "diff"]:
            return ""
        raise AssertionError(command)

    def build(self, name, fortran="-O3 -DNDEBUG", c="-O3 -DNDEBUG"):
        source = self.root / name
        build = source / "build"
        build.mkdir(parents=True)
        (source / "third_party/spral").mkdir(parents=True)
        (build / "CMakeCache.txt").write_text(
            f"CMAKE_HOME_DIRECTORY:INTERNAL={source}\n"
            "CMAKE_GENERATOR:INTERNAL=Unix Makefiles\n"
            "CMAKE_BUILD_TYPE:STRING=Release\n")
        for target, language, flags, artifact in [
                ("kls", "C", c, "libkls.a"),
                ("kls_spral_scaling", "Fortran", fortran, "libspral.a"),
                ("kls_bench", "C", c, "kls_bench")]:
            directory = build / f"CMakeFiles/{target}.dir"
            directory.mkdir(parents=True)
            (directory / "flags.make").write_text(
                f"# compile {language} with /compiler/{language}\n"
                f"{language}_FLAGS = {flags}\n"
                f"{language}_DEFINES = -DKLS_HAVE_SPRAL_SCALING=1\n"
                f"{language}_INCLUDES = -I{source}/include -I{build}/generated\n")
            obj = f"CMakeFiles/{target}.dir/source.o"
            (build / obj).write_bytes(b"object")
            command = (f"/compiler/C {obj} -o kls_bench libkls.a libspral.a -lm\n"
                       if target == "kls_bench" else f"/usr/bin/ar qc {artifact} {obj}\n")
            (directory / "link.txt").write_text(command)
            (build / artifact).write_bytes(b"artifact")
        return build / "kls_bench"

    def test_matching_builds_and_different_roots(self):
        proof = provenance.require_matching_builds(
            {"before": self.build("before"), "after": self.build("after")})
        self.assertEqual(proof["before"], proof["after"])
        self.assertEqual(len(proof["after"]["targets"]), 3)
        self.assertIn("Fortran", proof["after"]["compilers"])
        self.assertIn("<SOURCE>", str(proof))
        self.assertNotIn(str(self.root), str(proof))

    def test_c_matches_but_fortran_unoptimized(self):
        with self.assertRaisesRegex(ValueError, "Fortran_FLAGS"):
            provenance.require_matching_builds(
                {"before": self.build("before"), "after": self.build("after", fortran="")})

    def test_flag_order_is_not_normalized_away(self):
        with self.assertRaisesRegex(ValueError, "Fortran_FLAGS"):
            provenance.require_matching_builds({
                "before": self.build("before", fortran="-O0 -O3"),
                "after": self.build("after", fortran="-O3 -O0")})

    def test_missing_cache_fails_closed(self):
        before, after = self.build("before"), self.build("after")
        (after.parent / "CMakeCache.txt").unlink()
        with self.assertRaisesRegex(ValueError, "preflight failed"):
            provenance.require_matching_builds({"before": before, "after": after})

    def test_missing_dependency_flags_fails_closed(self):
        before, after = self.build("before"), self.build("after")
        (after.parent / "CMakeFiles/kls_spral_scaling.dir/flags.make").unlink()
        with self.assertRaisesRegex(ValueError, "kls_spral_scaling"):
            provenance.require_matching_builds({"before": before, "after": after})

    def test_stale_binary_fails_closed(self):
        before, after = self.build("before"), self.build("after")
        os.utime(after, ns=(1, 1))
        with self.assertRaisesRegex(ValueError, "stale"):
            provenance.require_matching_builds({"before": before, "after": after})

    def test_dependency_revision_mismatch_fails_closed(self):
        binary = self.build("before")
        original = self.command_output
        def changed(command, cwd=None):
            return "wrong" if command[:2] == ["git", "rev-parse"] else original(command, cwd)
        with patch.object(provenance, "output", side_effect=changed):
            with self.assertRaisesRegex(ValueError, "gitlink"):
                provenance.collect_build_provenance(binary)

    def test_reduce_only_does_not_require_build_metadata(self):
        with patch.object(sys, "argv", ["run_slim_regression.py", "--before", "missing",
                "--after", "missing", "--before-revision", "old", "--output", str(self.root),
                "--reduce-only"]), patch.object(runner, "require_matching_builds") as check, \
                patch.object(runner, "reduce_records", return_value={}) as reduce, \
                patch("builtins.print"):
            runner.main()
        check.assert_not_called()
        reduce.assert_called_once_with(self.root)


if __name__ == "__main__":
    unittest.main()
