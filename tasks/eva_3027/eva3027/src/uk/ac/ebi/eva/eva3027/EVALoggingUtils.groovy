package uk.ac.ebi.eva.eva3027

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.LoggerContext
import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.ConsoleAppender
import ch.qos.logback.core.FileAppender
import org.slf4j.Logger
import org.slf4j.LoggerFactory

class EVALoggingUtils {
    static Logger getLogger(Class aClass, boolean suppressConsoleOutput = false, String fileName = null) {
        LoggerContext lc = (LoggerContext) LoggerFactory.getILoggerFactory()
        PatternLayoutEncoder ple = new PatternLayoutEncoder()
        ple.setPattern("%date %level [%thread] %logger{10} [%file:%line] %msg%n")
        ple.setContext(lc)
        ple.start()
        def scriptLogger = (ch.qos.logback.classic.Logger)LoggerFactory.getLogger(aClass)

        if (!Objects.isNull(fileName)) {
            def fileAppender = new FileAppender()
            fileAppender.setFile(fileName)
            fileAppender.setAppend(true)
            fileAppender.setEncoder(ple)
            fileAppender.setContext(lc)
            fileAppender.start()
            scriptLogger.addAppender(fileAppender)
            scriptLogger.setAdditive(false)
        }
        // For interactive use on Emacs terminal, being able to output only to files is desirable
        // if voluminous output is expected that can likely be lost in terminal scrolling
        // and also likely to result in out-of-memory issues
        if (!suppressConsoleOutput) {
            def consoleAppender = new ConsoleAppender()
            consoleAppender.setEncoder(ple)
            consoleAppender.setContext(lc)
            consoleAppender.start()
            scriptLogger.addAppender(consoleAppender)
            scriptLogger.setAdditive(false)
        }

        scriptLogger.setLevel(Level.INFO)
        return scriptLogger
    }
}
