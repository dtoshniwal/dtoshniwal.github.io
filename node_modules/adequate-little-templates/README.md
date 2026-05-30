# Adequate Little Templates

[![npm version](https://img.shields.io/npm/v/adequate-little-templates.svg)](https://www.npmjs.com/package/adequate-little-templates)

A minimal, CSP-safe templating language. Unlike the micro templates I could find, this does not compile to a JavaScript function requiring `eval` or equivalent, so it will run in more contexts and be safer for untrusted input.

Currently: ~3KB gzipped.

## Quick Start

```bash
npm i adequate-little-templates
```

```javascript
import { render, compile, registerFunction } from 'adequate-little-templates';

// One-shot render
render('Hello {{ name }}!', { name: 'World' });
// -> "Hello World!"

// Compile for reuse
const template = compile('Hello {{ name }}!');
template({ name: 'World' }); // -> "Hello World!"
template({ name: 'Again' }); // -> "Hello Again!"
```

## Syntax

### Variable Interpolation

Output values, HTML escaped:

```handlebars
{{ title }}
{{ meta.author }}
{{ deeply.nested.value }}
```

### Raw Output

Output without HTML escaping (for trusted HTML like highlighted excerpts):

```
{{+ excerpt +}}
{{+ meta.description +}}
```

### Pipes and Functions

Transform values with pipes:

```handlebars
{{ title | lowercase }}
{{ title | uppercase }}
{{ title | truncate(50) }}
{{ title | truncate(50, "…") }}
{{ name | lowercase | truncate(20) }}
```

Or use function call syntax:

```handlebars
{{ truncate(title, 50) }}
{{ replace(title, "foo", "bar") }}
```

Pipes are sugar for function calls. These are equivalent:

```handlebars
{{ title | truncate(50) }}
{{ truncate(title, 50) }}
```

### Conditionals

```handlebars
{{#if author}}
  <span>By {{ author }}</span>
{{/if}}

{{#if image}}
  <img src="{{ image }}">
{{:else}}
  <div class="placeholder"></div>
{{/if}}

{{#if eq(status, "published")}}
  <span class="green">Published</span>
{{:else if eq(status, "draft")}}
  <span class="yellow">Draft</span>
{{:else}}
  <span class="gray">Unknown</span>
{{/if}}
```

### Iteration

```handlebars
{{#each items as item}}
  <li>{{ item.name }}</li>
{{/each}}

{{#each items as item, index}}
  <li>{{ index }}: {{ item.name }}</li>
{{/each}}

{{#each items as item}}
  <li>{{ item }}</li>
{{:else}}
  <li>No items found</li>
{{/each}}

{{#each items | limit(3) as item}}
  <li>{{ item.name }}</li>
{{/each}}
```

### Escaping Delimiters

Output literal `{{` with a backslash:

```handlebars
\{{ this is not a template tag }}
```

## Built-in Functions

### Comparison

| Function | Description |
|----------|-------------|
| `eq(a, b)` | Equal |
| `ne(a, b)` | Not equal |
| `gt(a, b)` | Greater than |
| `lt(a, b)` | Less than |
| `gte(a, b)` | Greater than or equal |
| `lte(a, b)` | Less than or equal |

### Boolean

| Function | Description |
|----------|-------------|
| `and(a, b, ...)` | Logical AND |
| `or(a, b, ...)` | Logical OR |
| `not(a)` | Logical NOT |

### String

| Function | Description |
|----------|-------------|
| `lowercase(s)` | Convert to lowercase |
| `uppercase(s)` | Convert to uppercase |
| `trim(s)` | Remove leading/trailing whitespace |
| `truncate(s, n)` | Truncate to n chars with "..." |
| `truncate(s, n, suffix)` | Truncate with custom suffix |
| `replace(s, find, rep)` | Replace all occurrences |

### Array

| Function | Description |
|----------|-------------|
| `limit(arr, n)` | Take first n items |
| `first(arr)` | First item |
| `last(arr)` | Last item |
| `length(arr)` | Array (or string) length |
| `join(arr, sep)` | Join with separator |

### Misc

| Function | Description |
|----------|-------------|
| `default(val, fallback)` | Use fallback if val is falsy |
| `safeUrl(url)` | Sanitize URL, blocking dangerous protocols |

### URL Sanitization

Use `safeUrl` to protect against `javascript:` and other dangerous URL schemes when rendering user-provided URLs:

```handlebars
<a href="{{ meta.url | safeUrl }}">{{ title }}</a>
```

Allowed:
- Relative URLs: `/path`, `./path`, `../path`, `?query`, `#hash`
- HTTP(S): `https://example.com`
- FTP: `ftp://example.com/file.txt`
- Email: `mailto:user@example.com`
- Phone: `tel:+1234567890`

Blocked (returns empty string):
- `javascript:alert('xss')`
- `data:text/html,...`
- `vbscript:...`
- Any other protocol

## Custom Functions

Register your own functions:

```javascript
import { registerFunction, render } from 'adequate-little-templates';

registerFunction('double', (n) => n * 2);
registerFunction('formatDate', (ts) => new Date(ts).toLocaleDateString());
registerFunction('pluralize', (count, singular, plural) =>
  count === 1 ? singular : plural
);

render('{{ double(5) }}', {});
// -> "10"

render('{{ count }} {{ pluralize(count, "item", "items") }}', { count: 3 });
// -> "3 items"

// Works with pipes too
render('{{ num | double }}', { num: 21 });
// -> "42"
```

## Truthiness

The following values are **falsy**:

- `null`
- `undefined`
- `false`
- `0`
- `""` (empty string)
- `NaN`
- `[]` (empty array)
- `{}` (empty object)

Everything else is **truthy**.

## Error Handling

Errors render inline for easy debugging:

```text
[Error: use #each for arrays]
[Error: cannot render object]
[Error: unknown foo()]
[Error: #each needs array]
[Error: join() needs 2 args]
```

## API

### Functions

```typescript
render(template: string, data: TemplateData): string
```
Parse and render a template in one step.

```typescript
compile(template: string): Template
```
Parse a template and return a reusable render function. Use this when rendering the same template multiple times.

```typescript
registerFunction(name: string, fn: CustomFn): void
```
Register a custom function for use in templates.

### TypeScript

If you need explicit type annotations, all types are exported:

```typescript
import { render, TemplateData, Value, Template, CustomFn } from 'adequate-little-templates';

const data: TemplateData = { title: 'Hello', count: 42 };
const myFn: CustomFn = (a, b) => Number(a) + Number(b);
```
