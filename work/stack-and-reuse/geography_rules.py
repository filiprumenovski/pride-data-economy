"""Country rules preserved from the original geography compiler."""
from pathlib import Path
from collections import Counter
import json, re, unicodedata, hashlib, html

def js(s,default=None):
    return json.loads(s) if isinstance(s,str) and s else ([] if default is None else default)
def norm(s):
    return re.sub(r'[^a-z0-9]','',unicodedata.normalize('NFKD',s or '').encode('ascii','ignore').decode().lower())
def clean(s):
    return html.unescape(re.sub('<[^>]+>','',s or '')).strip()
def digest(s):return hashlib.sha256(s.encode()).hexdigest()[:20] if s else None
iso=js((Path(__file__).resolve().parent/'reference/iso3166-1.json').read_text())['3166-1']
country_names={}
display={'United States of America':'United States','Korea, Republic of':'South Korea',
         'Russian Federation':'Russia','Taiwan, Province of China':'Taiwan',
         'Iran, Islamic Republic of':'Iran','Türkiye':'Turkey','Viet Nam':'Vietnam',
         'Bolivia, Plurinational State of':'Bolivia','Venezuela, Bolivarian Republic of':'Venezuela'}
for x in iso:
    name=display.get(x['name'],x['name'])
    country_names[norm(name)]=name
    for field in ['name','alpha_2','alpha_3','official_name','common_name']:
        if x.get(field):country_names[norm(x[field])]=name
aliases={'UK':'United Kingdom','England':'United Kingdom','Scotland':'United Kingdom',
 'U.S.A.':'United States','The United States':'United States','the United States of America':'United States',
 'United States of America (USA)':'United States','The Netherlands':'Netherlands',
 'Czech Republic':'Czechia','Brasil':'Brazil','Italia':'Italy','Deutschland':'Germany',
 'Suisse':'Switzerland','Schweiz':'Switzerland','South Korea':'South Korea',
 'Republic of Korea':'South Korea','Rep. of Korea':'South Korea','Hong Kong, China':'Hong Kong',
 'Hong Kong S.A.R., China':'Hong Kong','Macau':'Macao','Macau SAR China':'Macao',
 'Macedonia':'North Macedonia','Turkiye':'Turkey','Republic of Ireland':'Ireland',
 'Republic of China (Taiwan)':'Taiwan','Magyarorszag':'Hungary','Belgique':'Belgium',
 'Estados Unidos':'United States','Allemagne':'Germany','Nederland':'Netherlands',
 'Spagna':'Spain','Maroc':'Morocco','Tunisie':'Tunisia'}
country_names.update({norm(k):v for k,v in aliases.items()})
unknown_country=Counter()
def country(s):
    key=norm(s)
    if key in ('','na','null','none','unknown'):return None
    result=country_names.get(key)
    if result is None:unknown_country[s]+=1
    return result
assert country('USA')==country('U.S.')==country('United States')=='United States'
assert country('Korea') is None and country('CCMS') is None
assert country('Spagna')=='Spain' and norm('Université de Genève')==norm('Universite de Geneve')
address_names={display.get(x['name'],x['name']) for x in iso}|set(aliases)|{'USA','U.S.','U.S.A'}
address_names.update(x['official_name'] for x in iso if x.get('official_name'))
address_names.update(x['common_name'] for x in iso if x.get('common_name'))
address_names={s for s in address_names if len(s)>=4 or s in ['UK','US','USA','U.S.']}
address_names.discard('Georgia')
address_pattern=re.compile(r'(?:^|[,;\s])('+ '|'.join(re.escape(s) for s in sorted(address_names,key=len,reverse=True))+r')\s*[.)]?$',re.I)
def address_country(aff):
    aff=re.sub(r'\s*\(\s*lab head\s*\)\s*$','',aff or '',flags=re.I)
    match=address_pattern.search(aff)
    return country(match.group(1)) if match else None
assert address_country('BioArCh, University of York, UK')=='United Kingdom'
assert address_country('University of Michigan, USA')=='United States'
assert address_country('University of Michigan, United States of America')=='United States'
assert address_country('University of Georgia') is None

