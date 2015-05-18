# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import sys
import os
import random
import time
import re
import datetime
import cPickle as pickle
import argparse

from subrunner import tresult

user = os.environ.get("USER", "NOUSER")
tmppath = "/tmp/{0}_".format(user)
errfile = tmppath + "test_error"
resfile = tmppath + "test_result"

def absPath(relpath, packagename=""):
    """Dodgy but useful way to find an absolute path """
    
    if os.path.exists(relpath) or relpath[0]=="/": 
        return relpath
    
    path = os.path.abspath(sys.argv[0])
    return "{0}/{1}".format(path[0 : path.rfind(packagename) + len(packagename)], relpath)

class TestRunner(object):

    def __init__(self, testdirs=[], outfile=None, pattern=None,
                 test_list_file=None, seed=0, randomize=False):
        """
        :param testdirs: the directories in which to look for tests
        @type testdirs : list
        :param pattern: only runs tests containing pattern; e.g. if pattern="scrape" then only tests like test_scraper.py, test_scraperspeed.py will be run
        @type pattern : str
        :param test_list: run only the tests listed in the file given by test_list. Names of tests much each be on a separate line.
        @type test_list : str
        :param randomize: if True, run tests in a random order
        @type randomize : bool
        :param seed: the seed for the PRNG -- if a particular seed causes test failures, then it may be useful to pass that seed back in later for debugging.
        @type seed : int
        
        NB: both pattern and test_list are applied, so if pattern="cache" and 
        test_list contains the lines:
          test_cache
          test_excel
          test_normalization
        only test_cache will be run.
        """
        
        self.testdirs = testdirs
        self.pattern = pattern
        self.test_list = []
        self._readTestList(test_list_file)
        self.overwrite_sys = False
        self.randomize = randomize
        
        if not seed:
            self.seed = long(time.time())
        else:
            self.seed = seed
            
        random.seed(self.seed)
        
        if outfile:
            self.out = open(outfile, "w")
            self.saved_stdout = sys.stdout
            self.saved_stderr = sys.stdout
            sys.stdout = self.out
            sys.stderr = self.out
            self.overwrite_sys = True
        else: 
            self.out = sys.stderr

    def _readTestList(self, test_list_file):
        """Read a list of modules to test from a file containing one test
        module file name per line.

        :param test_list_file: name of file containing the test module list
        @type test_list_file : str
        """
        self.test_list = []
        
        try:
            with open(test_list_file, "r") as tl_handle:
                lines = tl_handle.readlines()
        except:
            return
        
        for l in lines:
            l = l.strip()
            l = l.replace(".py", "")
            self.test_list.append(l)

    def __del__(self):
        """Clean up on destruction: restore stdio to its prior state.
        """
        
        if self.overwrite_sys:
            sys.stdout = self.saved_stdout
            sys.stderr = self.saved_stderr
    
    def getTestModules(self):
        """ Locates all files in the testdirs of packagename of the form
        test_*.py, and imports the corresponding modules into a list which
        it iterates through, calling runTests() on each in turn while
        recording the results.
        """

        import tests
        allmodules = []
        
        for testdir in self.testdirs:
            files = os.listdir(absPath(testdir))
            test_re = re.compile("test_\w+\.py$", re.IGNORECASE)
            
            testfiles = [f for f in files if test_re.match(f)]
            
            testfiles.sort(reverse=True)
            filenameToModuleName = lambda f: os.path.splitext(f)[0]
            
            shortnames = map(filenameToModuleName, testfiles)
            modulenames = [(testdir + mod).replace("/", ".") for mod in shortnames]
            allmodules += modulenames
            
        #exclude all tests that don't contain pattern (e.g. "provenance") if
        #pattern is set
        if self.pattern: 
            allmodules = [m for m in allmodules if self.pattern in m]

        # exclude all tests not in test_list if test_list is set
        if self.test_list:
            allmodules = [m for m in allmodules if m in self.test_list]
            
        return allmodules
    
    def runTestModules(self, modulenames):
        """Execute the test in a completely separate subrunner process so
        that there can be no cross-talk between tests.

        :param modulenames: modules to test
        @type modulenames : list
        @return: results, test errors
        @rtype : tuple
        """

        cmd = "{0} subrunner.py".format(sys.executable)
        cmd += " {0}"
        
        results = []
        test_errors = []
        count = 0
        
        for name in modulenames:
            time.sleep(0.5) #without a sleep, testrunner is hard to Ctrl-C
            count += 1
            self.out.write("\n----- Starting test: {0} -----\n".format(name))
            runcmd = cmd.format(name)
	    # print "Run Cmd %s" %runcmd
            os.system(runcmd)
            with open(errfile, "r") as err_handle:
                err = pickle.load(err_handle)
            
            if err:
                self.out.write("{0}\n".format(err))
                test_errors.append((name, err))

            with open(resfile, "r") as res_handle:
                res = pickle.load(res_handle)
            results.append((name, res))
            
            self.out.write("----- Ending test: {0} ({1}/{2})-----\n\n".format(name, count, len(modulenames)))
            
        return results, test_errors
    
    def runlast(self,testfiles, lasttest):
        """Reorder test list to run lasttest last.

        :param testfiles: all test file modules to be run
        @type testfiles : list
        :param lasttest: name of test module to run last
        @type lasttest : str
        """
        
        try:
            testfiles.remove(lasttest)
            testfiles.append(lasttest)
        except ValueError:
            pass
    
    def printresults(self, results, test_errors):
        """Show test results, including failures and errors.

        :param results: summary results for tests
        @type results : list
        :param test_errors: information about severely failed test modules
        @type test_errors : list
        @return: 0 for a perfect run, 1 if test errors or failures encountered
        @rtype : int
        """
        
        self.out.write("\nSummary:\n")
        
        def rep(tr):
            try:
                testsRun = tr.testsRun
                testsRunNumber = tr.testsRun
            except AttributeError:
                testsRun = "FAILED -- test result is None"
                testsRunNumber = 0
            try:
                errors = len(tr.errors)
            except AttributeError:
                errors = 1
            try:
                failures = len(tr.failures)
            except AttributeError:
                failures = 0
                
            return ("run={0} errors={1} failures={2}".format(testsRun, errors,
                                                             failures),
                    (testsRunNumber, errors, failures))
        
        run, errors, failures = 0, 0, 0
        if not results:
            self.out.write("NO tests were run.\n")
            
        for result in results:
            result_string, stats = rep(result[1])
            self.out.write("{0}, {1}\n".format(result_string, result[0]))
            run += stats[0]
            errors += stats[1]
            failures += stats[2]
            
        if test_errors:
            self.out.write("\nThe following tests did NOT run:\n")
            for t in test_errors:
                self.out.write("{0} [{1}]\n".format((t[0],str(t[1]))))
            self.out.write("\n")
        self.out.write("Total run={0} errors={1} failures={2}, did not run={3}\n".format(run, errors, failures, len(test_errors)))
        
        if errors or failures or test_errors:
            self.out.write("FAILED\n")
            return 1
        
        else:
            return 0
    
    def run(self):
        modules = self.getTestModules()
        if self.randomize:
            random.shuffle(modules)
            
        mlist = ", ".join(modules)
        
        self.out.write("UNIT TEST RUN: {0}\n".format(datetime.datetime.now()))
        
        if self.randomize:
            self.out.write("Running tests in random order using seed: {0}.\n".format(self.seed))
        else:
            self.out.write("Running tests in reverse alphabetial order.\n")
            
        self.out.write("Running tests in the following modules: {0}\n".format((str(mlist))))
        results, test_errors = self.runTestModules(modules)
        
        return self.printresults(results, test_errors)
   
def othertr():
    """Check if this file's name appears in the process list. If so, another
    testrunner may be simultaneously active -- warn the user.

    @return: True if another testrunner may be in use
    @rtype : bool
    """
    
    procs = os.popen("ps -ef | grep {0}".format(__file__)).read() 
    procs = procs.split("\n")
    procs = [p for p in procs if p and "grep" not in p and "emacs" not in p]
    
    return len(procs) > 1

def usage():
    msg="""Usage: python {0} -o outfile -p pattern -t testlistfile -r randomize -s seed\n\nAll parameters are optional.\n\n\toutfile: the file to which all output is written; stderr is the default.\n\tpattern: if present only tests containing the substring pattern will be run \n\t\t\n\ttestlistfile: a file containing a list of specific tests to run, one per line \n\t\t(nb: both pattern and testlistfile can be used together).\n\trandomize: if 1, run the tests in a random order.\n\tseed: the seed for the randomizer.\n\n""".format(__file__)
    sys.stderr.write(msg)

if __name__ == "__main__":
    if othertr():
        sys.stderr.write("There appears to be another testrunner running.\n")
        time.sleep(3)
        
    try:
        if sys.argv[1] == "-?":
            usage()
            sys.exit(0)
    except IndexError:
        pass
    
    #parser = OptionParser()
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile", type=str, help="the file to which all output is written; stderr is the default")
    parser.add_argument("-p", "--pattern", type=str, help="if present only tests containing the given substring will be run")
    parser.add_argument("-t", "--test_list", type=str, help="name of a file containing a list of specific tests to run, one per line")
    parser.add_argument("-r", "--randomize", type=int, default=0, help="randomize order of test module execution, if value is not zero")
    
    options = parser.parse_args()
    # since the test files look at sys.argv, we need to clear it
    sys.argv = [sys.argv[0]]
    
    tr=TestRunner(outfile=options.outfile, pattern=options.pattern,
                  test_list_file=options.test_list,
                  randomize=options.randomize, seed=options.randomize,
                  testdirs=["tests/"])
    
    result = tr.run()
    sys.exit(result)

