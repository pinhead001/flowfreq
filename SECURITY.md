# Security Policy

## Reporting a Vulnerability

If you discover a security vulnerability in FlowFreq, **please do NOT open a public GitHub issue**.

Instead, please report it by email to **pinhead001@github.com** with the following information:

1. **Description** — What is the vulnerability and how can it be exploited?
2. **Affected versions** — Which releases of FlowFreq are impacted?
3. **Steps to reproduce** — Minimal code example (if applicable)
4. **Suggested fix** — (Optional) Your proposed solution
5. **CVSS score** — (Optional) Severity assessment

### What happens next

- We will acknowledge receipt **within 48 hours**
- We will provide an estimated timeline for a fix
- We will request a **90-day grace period** before public disclosure to allow users time to upgrade
- You will be credited in the security advisory (unless you request anonymity)

## Scope

✅ **In scope:**
- Vulnerabilities in FlowFreq source code (flowfreq/ directory)
- Vulnerabilities in vendored dependencies (vendor/peakfqr)
- Remote code execution, SQL injection, or similar attack vectors
- Sensitive data leakage

❌ **Out of scope:**
- USGS NWIS API downtime or changes
- Client configuration errors
- Network connectivity issues
- Denial of service via network flooding

## Security Best Practices

When using FlowFreq:

- Keep Python and dependencies up to date: `pip install --upgrade flowfreq`
- Validate input data, especially USGS gage IDs from untrusted sources
- Use regional skew parameters from authoritative USGS sources
- Do not commit API keys or credentials in code; use environment variables

## Vulnerability History

No security vulnerabilities have been reported as of 2026-08-31.
