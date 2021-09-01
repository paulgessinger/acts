#!/usr/bin/env python3
import json
import os
import sys

import asyncio
import aiohttp
from gidgethub.aiohttp import GitHubAPI


async def main():
    context = json.loads(os.environ["GITHUB_CONTEXT"])
    token = context["token"]
    repo = context["repository"]
    async with aiohttp.ClientSession(loop=asyncio.get_event_loop()) as session:
        gh = GitHubAPI(session, __name__, oauth_token=token)

        # get up to date pr info
        pr = await gh.getitem(
            f"repos/{repo}/pulls/{context['event']['pull_request']['number']}"
        )

        if pr["milestone"] is None:
            print("::error::No milestone is set on the PR")
            sys.exit(1)
        else:
            print("Milestone", pr["milestone"]["title"], "is set")
            sys.exit(0)


if __name__ == "__main__":
    loop = asyncio.get_event_loop()
    loop.run_until_complete(main())
