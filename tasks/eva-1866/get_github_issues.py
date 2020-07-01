import getpass

import sys
from cached_property import cached_property
from github import Github
from operator import itemgetter


class AllIssues:

    def __init__(self, github, user_name):
        self.g = github
        self.user_name = user_name
        self.milestones = None

    @cached_property
    def user_account(self):
        return self.g.get_user(self.user_name)

    @cached_property
    def repos(self):
        return [repo for repo in self.user_account.get_repos()]

    def _issues_per_repo(self, pr=False):
        issues_per_repo = {}
        for repo in self.repos:
            issues_per_repo[repo.name] = []
            for issue in repo.get_issues(state='open'):
                if bool(issue.pull_request) is bool(pr):
                    issues_per_repo[repo.name].append(issue)
        return issues_per_repo

    def _print_repo_issue_with_last_event(self, pr=False):
        issue_per_repo = self._issues_per_repo(pr)
        for repo in self.repos:
            issues = issue_per_repo.get(repo.name, pr)
            if not issues:
                sys.stderr.write('no issues for %s\n' % repo.name)
                continue

            sys.stderr.write('download %s issues for %s\n' % (len(issues), repo.name))
            for issue in issues:
                events = [(event.event, event.actor.login, event.created_at.strftime('%Y-%m-%d')) for event in issue.get_events()]
                events.sort(key=itemgetter(2))
                if events:
                    last_event_action, last_event_user, last_event_time = events[-1]
                else:
                    (last_event_action, last_event_user, last_event_time) = ('Creation', issue.user.login, issue.created_at.strftime('%Y-%m-%d'))
                print('\t'.join([
                    repo.name,
                    str(issue.number) + ': ' + issue.title,
                    issue.created_at.strftime('%Y-%m-%d'),
                    str(issue.user.login),
                    last_event_action,
                    last_event_user,
                    last_event_time,
                    issue.html_url
                ]))

    def print_repo_issue_with_last_event(self):
        self._print_repo_issue_with_last_event()

    def print_repo_pr_with_last_event(self):
        self._print_repo_issue_with_last_event(pr=True)


if __name__ == '__main__':
    password = getpass.getpass('Github password for username ' + getpass.getuser())
    g = Github(getpass.getuser(), password)
    all_issues = AllIssues(g, 'EBIvariation')
    all_issues.print_repo_issue_with_last_event()
    all_issues.print_repo_pr_with_last_event()
