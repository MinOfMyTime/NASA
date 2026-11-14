"""Probe Earthdata login and granule availability for the Earth products."""
import earthaccess

auth = earthaccess.login(strategy="netrc")
print("authenticated:", auth.authenticated)
for sn, v, t in [
    ("MYD11C3", "061", ("2023-01-01", "2023-12-31")),
    ("MCD43C3", "061", ("2023-01-01", "2023-12-31")),
    ("AG1KM", "003", None),
    ("AG100", "003", None),
    ("EMITL2BMIN", "001", None),
]:
    kw = dict(short_name=sn, version=v, count=3)
    if t:
        kw["temporal"] = t
    r = earthaccess.search_data(**kw)
    print(sn, len(r), r[0].data_links()[:2] if r else None, [round(g.size(), 2) for g in r])
