package eva3399

import org.slf4j.LoggerFactory

class Utils {
    static def logger = LoggerFactory.getLogger(Utils.class)

    static Object runProcess(String command) {
        logger.info("Running command: ${command}...")
        def process = ["bash", "-c", command].execute()
        def (out, err) = [new StringBuffer(), new StringBuffer()]
        process.waitForProcessOutput(out, err)
        if (!out.toString().trim().equals("")) logger.info(out.toString())
        if (process.exitValue() != 0) {
            logger.error("Command: $command exited with exit code ${process.exitValue()}!!")
            logger.error("See error messages below: \n" + err.toString())
            throw new Exception(err.toString())
        }
        logger.warn(err.toString())
        logger.info("Command: $command completed successfully!")
        return ["out": out.toString(), "err": err.toString()]
    }

    static String[] getCustomFastaAndAssemblyReportPaths (String assembly, String fastaDir) {
        def customFastaPath = runProcess("find ${fastaDir} -maxdepth 3 -iname ".toString() +
                "'${assembly}_custom.fa'| head -1".toString()).out.trim()
        def customAssemblyReportPath = runProcess("find ${fastaDir} -maxdepth 3 -iname ".toString() +
                "'${assembly}_assembly_report_custom.txt'| head -1".toString()).out.trim()
        if (customAssemblyReportPath.equals("")) throw new Exception("Could not find assembly report for assembly $assembly!!")
        return [customFastaPath, customAssemblyReportPath]
    }
}
