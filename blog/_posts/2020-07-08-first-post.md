---
layout: post
title: '"Sometimes the best way to learn something...'
---

... is by doing it wrong and then looking at what you did." - Neil Gaiman

Earlier this year I got introduced to [Observable](http://observablehq.com/), a D3-based platform for making interactive data visualizations. A few weeks ago I started making an interative spline editor there and I have to say that it has been an enjoyable learning experience.

I ended up filling up the notebook so that it can also serve as a gentle, gentle introduction to interactive design using smooth splines. If you have the time, [check it out here](https://observablehq.com/@dtoshniwal/interactive-design-with-smooth-splines).

<div id="observablehq-94ef5caf"></div>
<script type="module">
import {Runtime, Inspector} from "https://cdn.jsdelivr.net/npm/@observablehq/runtime@4/dist/runtime.js";
import define from "https://api.observablehq.com/@dtoshniwal/interactive-design-with-smooth-splines.js?v=3";
const inspect = Inspector.into("#observablehq-94ef5caf");
(new Runtime).module(define, name => name === "graphic" ? inspect() : undefined);
</script>

All in all, happy with what the platform can do; time permitting, future posts will diver deeper into splines or their applications!
