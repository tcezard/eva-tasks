import subprocess
import sys
from os.path import join
from time import sleep

from ebi_eva_common_pyutils.logger import AppLogger


class ScriptWriter(AppLogger):
    """
    Writes a basic job submission script. Subclassed by SlurmWriter. Initialises with self.lines as an empty list, which
    is appended by self.write_line. This list is then saved line by line to self.script_file by self.save.
    """
    suffix = '.sh'
    array_index = 'JOB_INDEX'
    mapping = {
        'cpus': '# cpus: {}',
        'exclusive': '# exclusive',
        'job_name': '# job name: {}',
        'job_queue': '# queue: {}',
        'jobs': '# job array: 1-{}',
        'log_file': '# log file: {}',
        'mem': '# mem: {}gb',
        'walltime': '# walltime: {}',
    }

    def __init__(self, job_name, working_dir, job_queue=None, log_commands=True, **kwargs):
        self.log_commands = log_commands
        self.working_dir = working_dir
        self.log_file = join(self.working_dir, job_name + '.log')
        self.parameters = dict(
            kwargs,
            job_name=job_name,
            job_queue=job_queue,
            log_file=self.log_file
        )

        self.script_name = join(working_dir, job_name + self.suffix)
        self.debug('Writing script: %s', self.script_name)
        self.lines = []

    def register_cmd(self, cmd, log_file=None):
        if log_file:
            cmd += ' > %s 2>&1' % log_file
        self.add_line(cmd)

    def register_cmds(self, *cmds, parallel):
        if parallel:
            self.add_job_array(*cmds)
        else:
            self.lines.extend(list(cmds))

    def add_job_array(self, *cmds):
        if self.parameters.get('jobs'):
            raise Exception('Already written a job array - can only have one per script')

        if len(cmds) == 1:
            self.register_cmd(cmds[0])
        else:
            self._start_array()
            for idx, cmd in enumerate(cmds):
                self._register_array_cmd(
                    idx + 1,
                    cmd,
                    log_file=self.log_file + str(idx + 1) if self.log_commands else None
                )
            self._finish_array()
            self.parameters['jobs'] = len(cmds)

    def _register_array_cmd(self, idx, cmd, log_file=None):
        """
        :param int idx: The index of the job, i.e. which number the job has in the array
        :param str cmd: The command to write
        """
        line = str(idx) + ') ' + cmd
        if log_file:
            line += ' > ' + log_file + ' 2>&1'
        line += '\n' + ';;'
        self.add_line(line)

    def add_line(self, line):
        self.lines.append(line)

    def _start_array(self):
        self.add_line('case $%s in' % self.array_index)

    def _finish_array(self):
        self.add_line('*) echo "Unexpected %s: $%s"' % (self.array_index, self.array_index))
        self.add_line('esac')

    def line_break(self):
        self.lines.append('')

    def save(self):
        """Save self.lines to self.script_name."""
        with open(self.script_name, 'w') as f:
            f.write('\n'.join(self.lines) + '\n')

    def add_header(self):
        """Write a header for a given resource manager. If multiple jobs, split them into a job array."""
        header_lines = ['#!/bin/bash\n']
        for k in sorted(self.parameters):
            v = self.parameters[k]
            header_lines.append(self.mapping[k].format(v))

        header_lines.extend(['', 'cd ' + self.working_dir, ''])

        # prepend the formatted header
        self.lines = header_lines + self.lines


class SlurmWriter(ScriptWriter):
    """Writes a Bash script runnable on Slurm"""
    suffix = '.slurm'
    array_index = 'SLURM_ARRAY_TASK_ID'
    mapping = {
        'cpus': '#SBATCH --cpus-per-task={}',
        'exclusive': '#SBATCH --exclusive',
        'job_name': '#SBATCH --job-name="{}"',
        'job_queue': '#SBATCH --partition={}',
        'jobs': '#SBATCH --array=1-{}',
        'log_file': '#SBATCH --output={}',
        'mem': '#SBATCH --mem={}g',
        'walltime': '#SBATCH --time={}:00:00'
    }


running_executors = {}


def stop_running_jobs():
    for job_id in running_executors:
        running_executors[job_id].cancel_job()

    for job_id in list(running_executors):
        running_executors[job_id].join()


class ClusterExecutor(AppLogger):
    script_writer = ScriptWriter
    finished_statuses = None
    unfinished_statuses = None
    submit_cmd = None

    def __init__(self, *cmds, prelim_cmds=None, pre_job_source=None, join_interval=None, **cluster_config):
        """
        :param cmds: Full path to a job submission script
        """
        self.interval = join_interval or 30
        self.job_id = None
        self.job_name = cluster_config.get('job_name')
        self.pre_job_source = pre_job_source
        self.cmds = cmds
        self.prelim_cmds = prelim_cmds
        self.writer = self.script_writer(**cluster_config)

    def write_script(self):
        if self.prelim_cmds:
            self.writer.register_cmds(*self.prelim_cmds, parallel=False)

        if self.pre_job_source:
            self.writer.register_cmd('source ' + self.pre_job_source)

        self.writer.line_break()
        self.writer.register_cmds(*self.cmds, parallel=True)
        self.writer.add_header()
        self.writer.save()

    def start(self):
        """Write the jobs into a script, submit it and capture qsub's output as self.job_id."""
        self.write_script()
        self._submit_job()
        running_executors[self.job_id] = self  # register to running_executors
        self.info('Submitted "%s" as job %s', self.writer.script_name, self.job_id)

    def join(self):
        """Wait until the job has finished, then return its exit status."""
        sleep(5)
        while not self._job_finished():
            sleep(self.interval)
        running_executors.pop(self.job_id, None)  # unregister from running_executors
        return self._job_exit_code()

    def _job_statuses(self):
        return ()

    def _job_exit_code(self):
        raise NotImplementedError

    def _submit_job(self):
        self.job_id = self._run_and_retry(self.submit_cmd + ' ' + self.writer.script_name)
        if self.job_id is None:
            raise Exception('Job submission failed')

    def _job_finished(self):
        statuses = self._job_statuses()
        for s in statuses:
            if s in self.finished_statuses:
                pass
            elif s in self.unfinished_statuses:
                return False
            else:
                raise Exception('Bad job status: %s' % s)
        return True

    def _get_stdout(self, cmd):
        p = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        exit_status = p.wait()
        o, e = p.stdout.read(), p.stderr.read()
        self.debug('%s -> (%s, %s, %s)', cmd, exit_status, o, e)
        if exit_status:
            return None
        else:
            return o.decode('utf-8').strip()

    def _run_and_retry(self, cmd, retry=3):
        attempt = 0
        while attempt < retry:
            msg = self._get_stdout(cmd)
            if msg is not None:
                return msg
            sleep(5)
            attempt += 1

    def cancel_job(self):
        if not self._job_finished():
            self._cancel_job()

    def _cancel_job(self):
        raise NotImplementedError


class SlurmExecutor(ClusterExecutor):
    script_writer = SlurmWriter
    unfinished_statuses = ('CONFIGURING', 'COMPLETING', 'PENDING', 'RUNNING', 'RESIZING', 'SUSPENDED',)
    finished_statuses = ('BOOT_FAIL', 'CANCELLED', 'COMPLETED', 'DEADLINE', 'FAILED', 'NODE_FAIL',
                         'PREEMPTED', 'TIMEOUT')
    submit_cmd = 'sbatch'

    def _submit_job(self):
        # sbatch stdout: "Submitted batch job {job_id}"
        super()._submit_job()
        self.job_id = self.job_id.split()[-1].strip()

    def _sacct(self, output_format):
        data = self._run_and_retry('sacct -nX -j {j} -o {o}'.format(j=self.job_id, o=output_format))
        return set(d.strip() for d in data.split('\n'))

    def _squeue(self):
        s = self._run_and_retry('squeue -h -j {j} -o %T'.format(j=self.job_id))
        if s:
            return set(s.split('\n'))

    def _job_statuses(self):
        s = self._squeue()
        if s:  # job is in squeue, so use that
            return s
        return set(s.rstrip('+') for s in self._sacct('State'))  # job no longer in squeue, so use sacct

    def _job_exit_code(self):
        exit_status = 0
        states = set()
        reports = self._sacct('State,ExitCode')
        for r in reports:
            state, exit_code = r.split()
            state = state.rstrip('+')
            states.add(state)
            exit_code = int(exit_code.split(':')[0])
            if state == 'CANCELLED' and not exit_code:  # cancelled jobs can still be exit status 0
                self.debug('Found a cancelled job - using exit status 9')
                exit_code = 9
            exit_status += exit_code

        self.info('Got %s states from %s (%s) with %s jobs: %s', len(states),
                  self.job_name, self.job_id, len(reports), states)
        return exit_status

    def _cancel_job(self):
        msg = self._run_and_retry('scancel ' + self.job_id)
        self.info(msg)


if __name__ == '__main__':
    commands = sys.stdin
    executor = SlurmExecutor(commands)
    executor.start()

