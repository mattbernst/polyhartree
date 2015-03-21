# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
import sys
import os
import cPickle as pickle

class tresult(object):
    def __init__(self): 
        self.failures = []
        self.errors = []
        self.testsRun = 0

def runtest():
    """Import a module into the local process for testing. Each test runs
    in a separate process to avoid side-effect contamination of other test
    modules in a suite.
    """
    
    err = ""
    result = ""
    user = os.environ.get("USER", "NOUSER")
    tmppath = "/tmp/{0}_".format(user)
    errfile = tmppath + "test_error"
    resfile = tmppath + "test_result"

    try:
        name = sys.argv[1]
        sys.argv = [sys.argv[0]]
        #trim off extra args that would alter test results by changing
        #test name or database loading behavior
        try:
            module = __import__(name, globals(), locals(), [""])
            res = module.runTests()
            sys.argv = [sys.argv[0]]
            
            #unless we change/remove these attributes, pickle will fail
            #because of explicitly or implicitly embdedded file objects
            del res.stream
            result = tresult()
            result.testsRun = res.testsRun
            
            for i in range(len(res.failures)):
                identifier, msg = res.failures[i]
                result.failures.append((str(identifier), msg))
                
            for i in range(len(res.errors)):
                identifier, msg = res.errors[i]
                result.errors.append((str(identifier), msg))
                
        except AttributeError:
            pass
        
        except Exception, e:
            msg = "There was an error {0} trying to run: {0}".format(e, name)
            err = msg
            
    finally:
        with open(errfile, "w") as errorhandle:
            pickle.dump(err, errorhandle)
        with open(resfile, "w") as resulthandle:
            pickle.dump(result, resulthandle)

if __name__ == "__main__":
    runtest()
