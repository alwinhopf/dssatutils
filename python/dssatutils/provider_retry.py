"""Bounded network backoff and bounded submission for large provider queues."""
import re
import time
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

class ProviderConnectivityError(RuntimeError):
    """Stop a provider stream without classifying points as no coverage."""

def provider_transient(error):
    return bool(re.search(r"resolve host|resolv.*host|name resolution|nodename nor servname|connection|connect to|timed? ?out|timeout|429|50[234]|temporarily unavailable", str(error), re.I))

def provider_retry(operation, attempts=3, delay_seconds=5, sleep=time.sleep):
    for attempt in range(attempts):
        try:
            return operation()
        except Exception as exc:
            if not provider_transient(exc):
                raise
            if attempt + 1 == attempts:
                raise ProviderConnectivityError(f"Provider connectivity unavailable after bounded retries: {exc}") from exc
            delay = min(delay_seconds * 2**attempt, 30)
            print(f"Provider connection failed; retry {attempt + 2}/{attempts} in {delay} seconds")
            sleep(delay)

def bounded_map(function, jobs, workers):
    """Submit only one job per worker ahead; cancel pending jobs on failure."""
    jobs = iter(jobs)
    if workers == 1:
        for job in jobs:
            yield function(job)
        return
    with ThreadPoolExecutor(max_workers=workers) as pool:
        pending = set()
        try:
            for _ in range(workers):
                job = next(jobs, None)
                if job is not None:
                    pending.add(pool.submit(function, job))
            while pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    yield future.result()
                    job = next(jobs, None)
                    if job is not None:
                        pending.add(pool.submit(function, job))
        finally:
            for future in pending:
                future.cancel()
