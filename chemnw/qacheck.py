# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import sys
import argparse
import os
import glob
import shlex
import subprocess
import stat
import pprint

class Verifier(object):
    def __init__(self):
        self.processed = []

    def generate_test_scripts(self, tests, directory,
                              top, target,
                              serial="runserial",
                              parallel="runmpi"):
        """Generate test battery scripts for serial and MPI execution from
        the given list of tests to be included. Each file is a C shell
        script similar to those bundled with NWChem QA but including only
        a subset of tests, ordered by increasing test cost.

        @param tests: tests to include in output, each entry (cost, name)
        @type tests : list
        @param directory: location of QA directory, where scripts will go
        @type directory : str
        @param top: NWCHEM_TOP environment variable
        @type top : str
        @param target: NWCHEM_TARGET environment variable
        @type target : str
        @param serial: name of serial-execution test battery script to generate
        @type serial : str
        @param parallel: name of MPI-execution test battery script to generate
        @type parallel : str
        """

        head = """#!/bin/csh -f
#
setenv NWCHEM_TOP {0}
setenv NWCHEM_TARGET {1}
setenv NWCHEM_EXECUTABLE `which nwchem`

set np = 1
if ($1 !="") then
  set np = $1
endif
""".format(top, target)

        sname = "{0}/{1}".format(directory, serial)
        pname = "{0}/{1}".format(directory, parallel)
        s = open(sname, 'w')
        p = open(pname, 'w')
        s.write(head)
        p.write(head)

        count = 0
        for cost, test in tests:
            count += 1
            tname = test.replace('.out', '')
            sline = "./runtests.unix procs $np {0} # estimated cost {1} number {2}\n".format(tname, cost, count)
            pline = sline.replace('runtests.unix', 'runtests.mpi.unix')
            s.write(sline)
            p.write(pline)

        s.close()
        p.close()

        for name in (sname, pname):
            st = os.stat(name)
            os.chmod(name, st.st_mode | stat.S_IEXEC)

    def find_ok_tests(self, qa_root, glob_pattern="doNightly*"):
        """Find test cases that already appear in QA test scripts
        included with NWChem. There are more test cases in
        the tests/ directory than those appearing in the bundled QA
        scripts, and those extra tests are often broken.

        N.B.: even many of the tests that do appear in bundled scripts
        are unreliable. Now only accepting tests that are in the doNightly
        script and that are not commented out.

        @param qa_root: path to the QA root directory
        @type qa_root : str
        @param glob_pattern: pattern to match test battery file names
        @type glob_pattern : str
        @return: names of test cases included in standard test batteries
        @rtype : set
        """

        cases = set()
        full_pattern = "{0}/{1}".format(os.path.abspath(qa_root), glob_pattern)
        for testfile in glob.glob(full_pattern):
            with open(testfile) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#'):
                        continue
                    if "runtests" in line:
                        tail = line.split("runtest", 1)[-1]

                        #This will get test case names plus some garbage
                        #like 's.unix'. The garbage doesn't matter because
                        #the names are only used to exclude bad tests.
                        for test in tail.split():
                            cases.add(test)

        return cases
                        

    def find_fast_tests(self, qa_root, core_seconds):
        """Find test cases that, according to the reference output files,
        executed in less than (n_cores * seconds) core_seconds. Some care
        must be taken to exclude broken tests. Tests are assumed not-broken
        only if they already appear in one of the QA do* scripts.

        @param qa_root: path to the QA root directory
        @type qa_root : str
        @param core_seconds: maximum number of core-seconds to include test
        @type core_seconds : int
        @return: tests that appear to run sufficiently fast
        @rtype : list
        """

        ok_tests = self.find_ok_tests(qa_root)

        tests = set()
        for root, dirs, files in os.walk(os.path.abspath(qa_root)):
            path = root.split(os.sep)
            #process only test reference outputs
            for f in files:
                if 'tests' in path and f.endswith('.out'):
                    fullfile = os.path.sep.join(path + [f])
                    nproc = 0
                    walltime = 0
                    with open(fullfile, 'r') as infile:
                        for line in infile:
                            line = line.strip()
                            if line.startswith('nproc'):
                                nproc = int(line.split()[-1])
                            elif line.startswith("Total times"):
                                s = line.split()[-1]
                                walltime = float(s[:-1])

                            if nproc and walltime:
                                cost = nproc * walltime
                                entry = (cost, f)
                                #last part of path is test directory,
                                #we're summing up costs of all tests
                                #in one directory
                                try:
                                    tests[path[-1]].add(entry)
                                except KeyError:
                                    tests[path[-1]] = set(entry)
                                    
        approved = []
        for k, v in tests.items():
            #skip tests that are too unreliable for standard QA scripts
            if k not in ok_tests:
                continue
            #test runner expects at least one file matching directory name
            files = set([x[1] for x in v])
            trial = "{0}.out".format(k)
            if trial not in files:
                continue

            cost = sum([x[0] for x in v])
            if cost < core_seconds:
                approved.append((cost, k))

        return approved
                    

    def parse_qa_log(self, logfile):
        """Extract tests run from a NWChem QA log file by fetching all
        'Running'/'verifying' pairs. Mark whether each test was OK or
        failed according to the log file. Failed tests can be examined
        more closely later.

        @param logfile: name of QA log file to open
        @type logfile : str
        """

        refs_path = os.path.dirname(os.path.abspath(logfile)) + '/testoutputs'

        with open(logfile, 'r') as f:
            for line in f.readlines():
                line = line.strip()

                if line.startswith('Running'):
                    test = line.split()[-1]
                    
                elif line.startswith('verifying'):
                    status = line.split()[-1]
                    test_name = test.split('/')[-1]
                    ref_file = "{0}/{1}.ok.out.nwparse".format(refs_path, test_name)
                    trial_file = "{0}/{1}.out.nwparse".format(refs_path, test_name)
                    entry = {'name' : test_name,
                             'basic_status' : status,
                             'reference' : ref_file,
                             'trial' : trial_file}
                    self.processed.append(entry)

    def execute(self, cmd, stdin_data=''):
        """Execute a command with subprocess.Popen, optionally supplying
        data to the command through stdin, and return the results.

        @param cmd: a command line for an external program
        @type cmd : str
        @param stdin_data: optional data to supply on stdin to external program
        @type stdin_data : str
        @return: data from stdout
        @rtype : str
        """

        command = shlex.split(cmd)
        with(open('/dev/null', 'w')) as devnull:
            p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                 stdin=subprocess.PIPE, stderr=devnull)
            output = p.communicate(input=stdin_data)[0]

        return output

    def numericize(self, line):
        """Split a line of text by whitespace and try to convert all
        numbers to floating point numeric values.

        e.g.

        'Root 1 singlet 19.48009 a.u. 530.08030 eV'

        becomes

        ['Root', 1.0, 'singlet', 19.48009, 'a.u.', 530.0803, 'eV']

        @param line: input line of text from a nwparse QA file
        @type line : str
        @return: mixed list of strings and floats
        @rtype : list
        """

        converted = []
        for piece in line.split():
            try:
                v = float(piece)
            except ValueError:
                v = piece
            converted.append(v)

        return converted

    def score_mismatch(self, reference, trial):
        """Score a trial output against a reference output. Reference and
        trial are both lines of text extracted from .nwparse files produced
        by a failing QA test.

        The final score is a tuple consisting of

        (gross_mismatch, numeric_absolute_difference)

        Higher scores indicate worse deviations. Gross mismatches are worse
        deviations than numerical differences.

        If lines match but the numeric values in them do not, the absolute
        difference of values is added to the numeric_absolute_difference.

        If lines do not even contain numeric values in corresponding positions,
        or one group has more lines than the other, the mismatches add to
        the gross_mismatch count.

        @param reference: mismatching lines from reference .nwparse file
        @type reference : list
        @param trial: mismatching lines from current QA trial .nwparse file
        @type trial : list
        @return: mismatch score, higher scores indicate worse mismatches
        @rtype : tuple
        """

        gross_mismatches = 0
        abs_numeric_diff = 0.0

        nlines = max(len(reference), len(trial))
        
        for k in range(nlines):
            try:
                r = self.numericize(reference[k])
                t = self.numericize(trial[k])

                ncolumns = max(len(r), len(t))
                
                for j in range(ncolumns):
                    try:
                        #if at least one line has a number in current column,
                        #try to get the numeric difference
                        if float in [type(r[j]), type(t[j])]:
                            diff = abs(r[j] - t[j])
                            abs_numeric_diff += diff
                            
                    #lines don't match in length or don't have numbers in
                    #the same position - gross mismatch
                    except (IndexError, TypeError):
                        gross_mismatches += 1
                        
            except IndexError:
                gross_mismatches += 1

        score = (gross_mismatches, abs_numeric_diff)
        return score

    def compare_outputs(self, reference_file, trial_file, attributes = None):
        """Parse reference and trial files and compare their contents.

        @param reference_file: a known-good NWChem log file for a calculation
        @type reference_file : str
        @param trial_file: an NWChem log file to compare against the reference
        @type reference_file : str
        """

        cmd = "diff {0} {1}".format(reference_file, trial_file)
        diff = self.execute(cmd).split('\n')

        # lines starting with < are differences in reference file,
        # lines starting with > are differences in trial file,
        # group differing lines and make them easy to compare
        r = [x[1:].strip() for x in diff if x.startswith('<')]
        t = [x[1:].strip() for x in diff if x.startswith('>')]
        
        if len(r) != len(t):
            #trial just has extra info not in ref, e.g. 'Total SCS-MP2 energy'
            if not r:
                mismatch_score = (0, 0.0)
                
            elif not t:
                base = os.path.basename(reference_file)
                import ipdb; ipdb.set_trace()

            else:
                mismatch_score = self.score_mismatch(r, t)

        else:
            mismatch_score = self.score_mismatch(r, t)

        return mismatch_score

    def compare_all_failures(self):
        """Perform a more detailed analysis of cases that did not pass the
        NWChem QA procedures.
        """

        failures = []
        for entry in self.processed:
            if entry['basic_status'] == 'failed':
                score = self.compare_outputs(entry['reference'],
                                             entry['trial'])
                if score != (0, 0.0):
                    failed = entry.copy()
                    failed['score'] = score
                    failures.append(failed)

        #sort from minor to major failures
        decorated = [(f['score'], f) for f in failures]
        decorated.sort()
        failures = [f[1] for f in decorated]
        
        pprint.pprint(failures)

        passed = len(self.processed) - len(failures)
        print ("Total {0} passed {1} failed {2}".format(len(self.processed), passed, len(failures)))


def main(args):
    v = Verifier()

    if args.test_root:
        approved = v.find_fast_tests(args.test_root, args.cost)
        approved.sort()
        v.generate_test_scripts(approved, os.path.abspath(args.test_root),
                                args.top, args.target)
        
    else:
        v.parse_qa_log(args.logfile)
        v.compare_all_failures()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Generate a NWChem test battery with serial and MPI execution scripts, or check the test results from a completed test battery.", epilog="Example: generate a fast test battery where each test has cost no greater than 100:\n qacheck.py -c 100 -t /opt/science/nwchem/Nwchem-6.3.revision25564-src.2014-05-03/QA\nExample: run tests and then check them: \n cd /opt/science/nwchem/Nwchem-6.3.revision25564-src.2014-05-03/QA\n ./runserial | tee quick.log\n qacheck.py -l quick.log")
    parser.add_argument('-c', '--cost', help="Maximum cost (wall clock time multiplied by number of processors) of tests to include in test battery.", type=int,default=1000)
    parser.add_argument('--top', help="NWCHEM_TOP location of tree where NWChem was built/installed.", default="/opt/science/nwchem/current")
    parser.add_argument('--target', help="NWCHEM_TARGET machine target that NWChem was built for.", default="LINUX64")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-t', '--test-root', help="Root directory location of NWChem QA, required to generate test battery, expressed either absolutely or relatively, e.g. '.' or '/opt/science/nwchem/Nwchem-6.3.revision25564-src.2014-05-03'.")
    group.add_argument('-l', '--logfile', help="NWChem QA log file from a test battery run, required to check QA output.")
    args = parser.parse_args()
    main(args)

