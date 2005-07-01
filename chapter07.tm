<TeXmacs|1.0.5>

<project|multigrid.tm>

<style|<tuple|article|skript|number-long-article>>

<\body>
  <assign|section-nr|6><section|Multiplicative Schwarz Methods>

  <\theorem>
    For the successive subspace correction method, we have

    <\equation*>
      <enorm|(id-P<rsub|0>)(id-P<rsub|1>)\<cdots\>(id-P<rsub|N>)||2>=1-<frac|1|1+C<rsub|0>>\<less\>1<space|1em>(C<rsub|0>\<neq\>0)
    </equation*>

    with

    <\equation*>
      C<rsub|0>=sup<rsub|<enorm|v||>=1>inf<rsub|v=<big|sum>w<rsub|n>><big|sum><rsub|n=0><rsup|N><enorm|P<rsub|n><big|sum><rsub|k=n+1><rsup|N>w<rsub|k>||2>
    </equation*>

    for <with|mode|math|w<rsub|n>\<in\>W<rsub|n>>, considering

    <\equation*>
      <enorm|P<rsub|n><big|sum><rsub|k=n+1><rsup|N>w<rsub|k>||2>=<enorm|P<rsub|0>(v-w<rsub|0>)||2>+<enorm|P<rsub|0>(v-w<rsub|0>-w<rsub|1>)||2>+<enorm|P<rsub|0>(v-w<rsub|0>-w<rsub|1>-w<rsub|2>)||2>.
    </equation*>
  </theorem>

  <\proof>
    <em|Step 1.> <with|mode|math|E<rsub|1>=id>.
    <with|mode|math|E<rsub|n>=E<rsub|n-1>(id-P<rsub|n>)>. Now,

    <\equation>
      <label|eq:muschwa-diff><enorm|v||2>-<enorm|E<rsub|N>v||2>=<big|sum><rsub|n=0><rsup|N>a(h*E<rsub|n-1>v,E<rsub|n-1>v)=<big|sum><rsub|n=0><rsup|N><enorm|P<rsub|n>E<rsub|n-1>v||2>
    </equation>

    and

    <\equation>
      <label|eq:muschwa-idsum><big|sum><rsub|m=0><rsup|n-1>P<rsub|m>E<rsub|m-1>+E<rsub|n-1>=id.
    </equation>

    We define

    <\eqnarray*>
      <tformat|<table|<row|<cell|<wide|I|~>:V>|<cell|\<rightarrow\>>|<cell|V<rsup|N+1>,<space|1em><wide|I|~>v\<assign\>(v,\<ldots\>,v)<rsup|T>>>|<row|<cell|<wide|E|~>:V>|<cell|\<rightarrow\>>|<cell|V<rsup|N+1>,<space|1em><wide|E|~>v\<assign\>(v,E<rsub|1>v,E<rsub|2>v,\<ldots\>,E<rsub|N-1>v)<rsup|T>>>|<row|<cell|<wide|L|~>:V<rsup|N+1>>|<cell|\<rightarrow\>>|<cell|V<rsup|N+1>,<space|1em><wide|L|~>\<assign\><matrix|<tformat|<table|<row|<cell|id>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|P<rsub|0>>|<cell|\<ddots\>>|<cell|>|<cell|>|<cell|>>|<row|<cell|\<vdots\>>|<cell|P<rsub|1>>|<cell|\<ddots\>>|<cell|>|<cell|>>|<row|<cell|\<vdots\>>|<cell|\<vdots\>>|<cell|\<ddots\>>|<cell|\<ddots\>>|<cell|>>|<row|<cell|P<rsub|0>>|<cell|P<rsub|1>>|<cell|\<cdots\>>|<cell|P<rsub|n-1>>|<cell|id>>>>>>>|<row|<cell|<wide|U|~>:V<rsup|N+1>>|<cell|\<rightarrow\>>|<cell|V<rsup|N+1>,<space|1em><wide|U|~>\<assign\><matrix|<tformat|<table|<row|<cell|0>|<cell|>|<cell|>|<cell|>>|<row|<cell|id>|<cell|\<ddots\>>|<cell|>|<cell|>>|<row|<cell|\<vdots\>>|<cell|\<ddots\>>|<cell|\<ddots\>>|<cell|>>|<row|<cell|id>|<cell|\<cdots\>>|<cell|id>|<cell|0>>>>>>>|<row|<cell|<wide|P|~>:V<rsup|N+1>>|<cell|\<rightarrow\>>|<cell|V<rsup|N+1>,<space|1em><wide|P|~>\<assign\>diag({P<rsub|n>}<rsub|n>)>>>>
    </eqnarray*>

    Note that <with|mode|math|<wide|L|~>> is invertible and
    <with|mode|math|<wide|U|~>> is nilpotent (such that, by the Neumann
    series, <with|mode|math|1\<pm\><wide|U|~>> is invertible). We observe
    <with|mode|math|<wide|L|~><wide|P|~>-<wide|P|~>=<wide|U|~><wide|P|~>>.
    This and (<reference|eq:muschwa-idsum>) is equivalent to
    <with|mode|math|<wide|L|~><wide|E|~>=<wide|I|~>>, and thus
    <with|mode|math|<wide|E|~>=<wide|L|~><rsup|-1><wide|I|~>>. Futhermore,
    (<reference|eq:muschwa-diff>) means that

    <\equation*>
      <enorm|v||2>-<enorm|E<rsub|N>v||2>=<enorm|<wide|P|~><wide|E|~>v|N+1|2>=<enorm|<wide|P|~><wide|L|~><rsup|-1><wide|I|~>v|N+1|2>
    </equation*>

    The norm is defined in an obvious manner:

    <\equation*>
      <enorm|(v<rsub|0>,\<ldots\>,v<rsub|T>)<rsup|T>|N+1|2>=<big|sum><enorm|v<rsub|n>||2>.
    </equation*>

    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|>|<cell|<enorm|(id-P<rsub|0>)\<cdots\>(id-P<rsub|N>)||2>>>|<row|<cell|>|<cell|=>|<cell|<enorm|E<rsub|N>||2>=sup<rsub|<enorm|v||>=1><enorm|E<rsub|N>v||2>=sup<rsub|<enorm|v||>=1><enorm|v||2>-<enorm|<wide|P|~><wide|L|~><rsup|-1><wide|I|~>v|N+1|2>>>|<row|<cell|>|<cell|=>|<cell|1-inf<rsub|<enorm|v||>=1><enorm|<wide|P|~><wide|L|~><rsup|-1><wide|I|~>v|N+1|2>.>>|<row|<cell|>|<cell|<above|=|(\<ast\>)>>|<cell|1-inf<rsub|<enorm|v||>=1><enorm|<wide|P|~><wide|L|~><rsup|-1><wide|P|~><wide|I|~>v|N+1|2>>>|<row|<cell|>|<cell|=>|<cell|1-inf<rsub|<enorm|v||>=1><enorm|<wide|P|~>(<wide|P|~>+<wide|V|~>)<rsup|-1><wide|P|~><wide|I|~>v|N+1|2>,>>>>
    </eqnarray*>

    where the step <with|mode|math|(\<ast\>)> uses

    <\eqnarray*>
      <tformat|<table|<row|<cell|<wide|L|~><rsup|-1>>|<cell|=>|<cell|(<wide|I|~>-(<wide|I|~>-<wide|L|~><rsup|>))<rsup|-1><below|=|<with|mode|text|Neumann>><big|sum><rsub|k=0><rsup|N>(<wide|I|~>-<wide|L|~>)<rsup|k>=<big|sum><rsub|k=0><rsup|N>(-<wide|U|~><wide|P|~>)<rsup|k><wide*|<wide|P|~>|\<wide-underbrace\>><rsub|<with|mode|text|can
      be added>>>>>>
    </eqnarray*>

    and we defined <with|mode|math|<wide|V|~>\<assign\><wide|P|~><wide|U|~><wide|P|~>>.

    Now, we show that

    <\equation*>
      <left|(>inf<rsub|<enorm|v||>=1><enorm|<wide|P|~>(<wide|P|~>+<wide|V|~>)<rsup|-1><wide|P|~><wide|I|~>v|N+1|2><right|)><rsup|-1>=1+C<rsub|0>,
    </equation*>

    with

    <\equation*>
      C<rsub|0>\<assign\>sup<rsub|<enorm|v||>=1><enorm|<wide|P|~><wide|U|~><rsup|T><wide|P|~><wide|I|~>v|N+1|>.
    </equation*>

    To start out, consider that <with|mode|math|(<wide|U|~>+<wide|U|~><rsup|T>+<wide|id|~>)=(N+1)<wide|id|~>>.
    Consider

    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|>|<cell|inf<rsub|v\<neq\>0><frac|N+1|<enorm|<wide|I|~>v|N+1|2>><enorm|<wide|P|~>(<wide|P|~>+<wide|V|~>)<rsup|-1><wide|P|~><wide|I|~>v|N+1|2>>>|<row|<cell|>|<cell|=>|<cell|inf<rsub|v\<neq\>0><frac|N+1|a(<wide|I|~>v,<wide|I|~>v)>a(<wide|P|~>(<wide|P|~>+<wide|V|~><rsup|T>)<rsup|-1><wide|P|~>(<wide|P|~>+<wide|V|~>)<rsup|-1><wide|P|~><wide|I|~>v,<wide|P|~><wide|I|~><wide|v|~>)>>|<row|<cell|>|<cell|=>|<cell|inf<rsub|v\<neq\>0>(N+1)<frac|a(<wide|P|~><wide|I|~>v,<wide|P|~><wide|I|~>v)|a(<wide|P|~>(<wide|P|~>+<wide|V|~>)<wide|P|~>(<wide|P|~>+<wide|V|~><rsup|T>)<wide|P|~><wide|I|~>v,<wide|P|~><wide|I|~>v)>>>|<row|<cell|>|<cell|=>|<cell|inf<rsub|v\<neq\>1>*<frac|(N+1)a(<wide|P|~><wide|I|~>v,<wide|P|~><wide|I|~>v)|(N+1)a(<wide|P|~><wide|I|~>v,<wide|P|~><wide|I|~>v)+a(<wide|P|~><wide|V|~><wide|V|~><rsup|T><wide|P|~><wide|I|~>v,<wide|P|~><wide|I|~>v)>>>|<row|<cell|>|<cell|=>|<cell|inf<rsub|(N+1)<enorm|<wide|P|~><wide|I|~>v|N+1|2>\<neq\>1><frac|1|1+<enorm|<wide|P|~><wide|V|~><rsup|T><wide|P|~><wide|I|~>v|N+1|2>>>>>>
    </eqnarray*>

    It remains to show that

    <\equation*>
      sup<rsub|(N+1)<enorm|<wide|P|~><wide|I|~>v|N+1|2>><enorm|<wide|P|~><wide|V|~><rsup|T><wide|P|~><wide|I|~>v|N+1|>=sup<rsub|<enorm|v||>=1>inf<rsub|v=<big|sum>w<rsub|n>><enorm|<wide|P|~><wide|U|~><rsup|T>w|N+1|2>
    </equation*>

    with <with|mode|math|w=(w<rsub|0>,\<ldots\>,w<rsub|N>)<rsup|T>>.
  </proof>

  <subsection|V-cycle estimate without regularity>

  Let <with|mode|math|V<rsub|0>\<subset\>V<rsub|1>\<subset\>\<cdots\>\<subset\>V<rsub|J>>
  be a nested sequence of finite element spaces

  <\eqnarray*>
    <tformat|<table|<row|<cell|a:V<rsub|J>\<times\>V<rsub|J>>|<cell|\<rightarrow\>>|<cell|\<bbb-R\>,<space|1em><with|mode|text|sym.pos.def.>>>|<row|<cell|<ip|\<cdot\>|\<cdot\>||>:V<rsub|J>\<times\>V<rsub|J>>|<cell|\<rightarrow\>>|<cell|\<bbb-R\>>>>>
  </eqnarray*>

  and

  <\eqnarray*>
    <tformat|<table|<row|<cell|<enorm|v||>>|<cell|\<assign\>>|<cell|<sqrt|a(v,v)>>>|<row|<cell|<l2norm|v||>>|<cell|\<assign\>>|<cell|<sqrt|<ip|v|v||>>>>>>
  </eqnarray*>

  with <with|mode|math|<ip|A<rsub|j>v<rsub|j>|w<rsub|j>||>=a(v<rsub|j>,w<rsub|j>)>,
  <with|mode|math|a(P<rsub|j>v,v<rsub|j>)=a(v,v<rsub|j>)>,
  <with|mode|math|<ip|Q<rsub|j>v|v<rsub|j>||>=<ip|v|v<rsub|j>||>>. We need a
  <em|stable decomposition>

  <\equation>
    <label|eq:stable-decomp><enorm|Q<rsub|0>v||2>+<big|sum>\<lambda\><rsub|j><l2norm|Q<rsub|j>v-Q<rsub|j-1>v||2>\<leqslant\>K<rsub|0><enorm|v||2>
  </equation>

  [Remember:

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|sum><l2norm|w<rsub|j>||2>>|<cell|\<leqslant\>>|<cell|C<l2norm|w||2>,>>|<row|<cell|<big|sum><enorm|w<rsub|j>||2>>|<cell|\<leqslant\>>|<cell|<wide|C|~><enorm|w||2>>>>>
  </eqnarray*>

  and

  <\equation*>
    v=<big|sum><rsub|j=0><rsup|J>(Q<rsub|j>-Q<rsub|j-1>)v,
  </equation*>

  using <with|mode|math|Q<rsub|J>=id>, <with|mode|math|Q<rsub|0>=0>.] and a
  ``<em|strengthened Cauchy-Schwarz inequality>''

  <\equation>
    <label|eq:strengthened-cs>a(v<rsub|j>,v<rsub|k>)\<leqslant\>K<rsub|1>\<varepsilon\><rsup|k-j><enorm|v<rsub|j>||>\<lambda\><rsub|k><rsup|1/2><l2norm|v<rsub|k>||>
  </equation>

  for <with|mode|math|k\<geqslant\>j>, <with|mode|math|v<rsub|j>\<in\>V<rsub|j>>,
  <with|mode|math|v<rsub|k>\<in\>V<rsub|k>> with a constant
  <with|mode|math|\<varepsilon\>\<in\>(0,1)>. With full regularity, both
  (<reference|eq:stable-decomp>) and (<reference|eq:strengthened-cs>) are
  easy. Using the smoothing theory from Chapter <with|color|red|XXXXX>6, we
  assume a sort-of decomposition

  <\equation*>
    V<rsub|j>=<big|sum><rsub|n=1><rsup|N<rsub|j>>W<rsub|n><rsup|j>,
  </equation*>

  and

  <\eqnarray*>
    <tformat|<table|<row|<cell|P<rsub|n><rsup|j>:V<rsub|j>>|<cell|\<rightarrow\>>|<cell|W<rsub|n><rsup|j>>>|<row|<cell|Q<rsub|n><rsup|j>:V<rsub|j>>|<cell|\<rightarrow\>>|<cell|W<rsub|n><rsup|j>>>>>
  </eqnarray*>

  with the smoother

  <\equation*>
    (id-R<rsub|j>A<rsub|j>)=<big|prod><rsub|n=1><rsup|N<rsub|j>>(id-P<rsub|n><rsup|j>).
  </equation*>

  Our new V-cycle algorithm is

  <\itemize>
    <item><with|mode|math|u<rsub|J>> is gives
    <with|mode|math|r<rsub|J>=b<rsub|J>-A<rsub|J>u<rsub|J>>

    <item>For <with|mode|math|j=J,\<ldots\>,1>:

    <\itemize>
      <item><with|mode|math|c<rsub|J=0>>

      <item>For <with|mode|math|n=1,\<ldots\>,N<rsub|j>>:

      <\itemize>
        <item><with|mode|math|c<rsub|j><rsup|n>=(A<rsub|n><rsup|j>)<rsup|-1>r<rsub|j>>

        <item><with|mode|math|c<rsub|j>\<assign\>c<rsub|j>+c<rsub|j><rsup|n>>

        <item><with|mode|math|r<rsub|j>\<assign\>r<rsub|j>-A<rsub|j>c<rsub|j>>
      </itemize>

      <item><with|mode|math|r<rsub|j-1>=Q<rsub|j-1>r<rsub|j>>
    </itemize>

    <item><with|mode|math|c<rsub|0>=A<rsup|-1>r<rsub|0>>

    <item><with|mode|math|u<rsub|J>\<assign\>u<rsub|J>+<big|sum><rsub|j>c<rsub|j>>.
  </itemize>

  While <with|mode|math|(id-B<rsub|j>A<rsub|j>)> remains undefined for now,
  we get

  <\equation*>
    id-B<rsub|J>A<rsub|J>=(id-P<rsub|0>)<big|prod><rsub|j=1><rsup|J><big|prod><rsub|n=1><rsup|N<rsub|j>>(id-P<rsub|n><rsup|j>).
  </equation*>

  Also, we demand an <with|mode|math|L<rsup|2>>-estimate of the approximation
  quality (cf. <inactive|<reference|lem:2.12>>)

  <\equation>
    <label|eq:l2-approx-quality>\<lambda\><rsub|j><l2norm|v<rsub|j>-Q<rsub|j-1>v<rsub|j>||2>\<leqslant\>K<rsub|2><enorm|v<rsub|j>||2>.
  </equation>

  Assume that we have a partition of unity
  <with|mode|math|<big|sum><rsub|n=1><rsup|N<rsub|j>>\<psi\><rsub|n><rsup|j>\<equiv\>1>
  such that

  <\equation*>
    W<rsub|n><rsup|j>={v\<in\>V<rsub|j>:supp v\<subset\><wide*|supp
    \<psi\><rsub|n><rsup|j>|\<wide-underbrace\>><rsub|<wide|\<Omega\>|\<bar\>><rsub|n><rsup|j>>}.
  </equation*>

  We claim that

  <\enumerate-alpha>
    <item><with|mode|math|<l2norm|\<nabla\>\<psi\><rsub|n><rsup|j>|\<infty\>|2>\<leqslant\>C<rsub|0>\<lambda\><rsub|j>>.

    <item><with|mode|math|<enorm|\<Pi\><rsub|j>(\<psi\><rsub|n><rsup|j>w<rsub|j>)||2>\<leqslant\>C<rsub|1><enorm|\<psi\><rsub|n><rsup|j>w<rsub|j>||2>>,
    where <with|mode|math|\<Pi\><rsub|j>v=<big|sum><rsub|z\<in\>\<cal-N\><rsub|j>>v(z)\<varphi\><rsub|z><rsup|j>>
    and <with|mode|math|V<rsub|j>=span{\<varphi\><rsub|z><rsup|j>}>, where
    the <with|mode|math|\<varphi\><rsub|z><rsup|j>> are the nodal functions.
    (<em|Watch out:> <with|mode|math|\<Pi\><rsub|j>> is not a product over
    <with|mode|math|j>. It's the interpolation operator for the
    <with|mode|math|j>th level.) In the model case, we have

    <\equation*>
      <enorm|\<Pi\><rsub|j>(\<psi\><rsub|n><rsup|j>w<rsub|j>)||2>\<leqslant\>C<rsub|1><left|(>C<rsub|2><l2norm|\<psi\><rsub|n><rsup|j>|\<infty\>|><enorm|w<rsub|j>||2>+C<rsub|3><l2norm|\<nabla\>\<psi\><rsub|n><rsup|j>|\<infty\>|><l2norm|w<rsub|j>||2><right|)>.
    </equation*>

    <item><with|mode|math|<big|sum><rsub|j=1><rsup|N<rsub|j>><l2norm|v<rsub|j>|\<Omega\><rsub|n><rsup|j>|2>\<leqslant\>C<rsub|4><l2norm|v<rsub|j>||2>>
    and <with|mode|math|<big|sum><rsub|j=1><rsup|N<rsub|j>><enorm|v<rsub|j>|\<Omega\><rsub|n><rsup|j>|2>\<leqslant\>C<rsub|5><enorm|v<rsub|j>||2>>.
  </enumerate-alpha>

  <\theorem>
    <\equation*>
      <enorm|id-B<rsub|J>A<rsub|J>||2>=1-<frac|1|1+c<rsub|0>>,
    </equation*>

    where

    <\equation*>
      c<rsub|0>\<assign\>sup<rsub|<enorm|v||>=1>inf<rsub|v=v<rsub|0>+<big|sum><rsub|j=1><rsup|J><big|sum><rsub|n=1><rsup|N<rsup|j>>w<rsub|n><rsup|j>><enorm|P<rsub|0>(v-v<rsub|0>)||2>+<big|sum><rsub|j=1><rsup|J><big|sum><rsub|n=1><rsup|N<rsup|j>><enorm|P<rsub|n><rsup|j><big|sum><rsub|(k,m)\<gtr\>(j,n)>w<rsub|n><rsup|k>||2>.
    </equation*>
  </theorem>

  <\proof>
    <with|mode|math|v<rsub|0>=Q<rsub|0>v>,
    <with|mode|math|v<rsub|j>=(Q<rsub|j>-Q<rsub|j-1>)v>,
    <with|mode|math|v<rsub|j>=<big|sum>w<rsub|n><rsup|j>> with
    <with|mode|math|w<rsub|n><rsup|j>=\<Pi\><rsub|j>(\<psi\><rsub|n><rsup|j>v<rsub|j>)>,
    so that

    <\equation*>
      v<rsub|j>=<big|sum>w<rsub|n><rsup|j>=\<Pi\><rsub|j><left|(><big|sum><rsub|n=1><rsup|N<rsub|j>>\<psi\><rsub|n><rsup|j>v<rsub|j><right|)>.
    </equation*>

    <\eqnarray*>
      <tformat|<table|<row|<cell|c<rsub|0>>|<cell|\<leqslant\>>|<cell|<wide*|<enorm|P<rsub|0>(v-v<rsub|0>)||2>|\<wide-underbrace\>><rsub|<sqrt|\<cdot\>>\<leqslant\><enorm|v-v<rsub|0>||>\<leqslant\><enorm|v||>+<enorm|Q*v<rsub|0>||>>+<big|sum><rsub|j=1><rsup|J><big|sum><rsub|n=1><rsup|N<rsub|j>><enorm|P<rsub|n><rsup|j><left|(><big|sum><rsub|m=n+1><rsup|N<rsub|j>>w<rsub|m><rsup|j>+<big|sum><rsub|k=j+1><rsup|J>v<rsub|j><right|)>||2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|(<enorm|v||>+<enorm|Q*v<rsub|0>||>)<rsup|2>+<big|sum><rsub|j=1><rsup|J><big|sum><rsub|n=1><rsup|N<rsub|j>><left|(><enorm|P<rsub|n><rsup|j><big|sum><rsub|m=n+1><rsup|N<rsub|j>>w<rsub|m><rsup|j>||>+<enorm|P<rsub|n><rsup|j><big|sum><rsub|k=j+1><rsup|J>v<rsub|j>||><right|)><rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|(<enorm|v||>+<enorm|Q*v<rsub|0>||>)<rsup|2>+<big|sum><rsub|j=1><rsup|J><big|sum><rsub|n=1><rsup|N<rsub|j>><left|(><enorm|<big|sum><rsub|m=m+1><rsup|N<rsub|j>>\<Pi\><rsub|j>(\<psi\><rsub|n><rsup|j>v<rsub|j>)||>+<wide|C|~><rsub|0><enorm|v<rsub|j>|\<Omega\><rsub|n><rsup|j>|2>+<wide|C|~><rsub|1>\<lambda\><rsub|j><enorm|v<rsub|j>|\<Omega\><rsub|n><rsup|j>|2><right|)><rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|(<enorm|v||>+<enorm|Q*v<rsub|0>||>)<rsup|2>+<big|sum><rsub|j=1><rsup|J><big|sum><rsub|n=1><rsup|N<rsub|j>><left|(><enorm|<big|sum><rsub|m=m+1><rsup|N<rsub|j>>\<Pi\><rsub|j>(\<psi\><rsub|n><rsup|j>v<rsub|j>)||>+<wide|C|~><rsub|2><l2norm|v<rsub|j>|\<Omega\><rsub|n><rsup|j>|><right|)><rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|2<enorm|v||2>+2<enorm|Q<rsub|0>v<rsub|0>||>+<big|sum><rsub|j=1><rsup|J>C<rsub|5><wide|C|~><rsub|0>\<lambda\><rsub|j><l2norm|Q<rsub|j>v<rsub|j>-Q<rsub|j-1>v<rsub|j>||2>+<big|sum><rsub|j=1><rsup|J><big|sum><rsub|n=1><rsup|N<rsub|j>><enorm|P<rsub|n><rsup|j>(v-Q<rsub|j>v)||2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<wide|C|~><rsub|3>K<rsub|0><enorm|v||2>+<big|sum><rsub|j=1><rsup|J><big|sum><rsub|n=1><rsup|N<rsub|j>><enorm|P<rsub|n><rsup|j>(v-Q<rsub|j>v)||2>>>>>
    </eqnarray*>

    <with|color|red|Quadrate aufräumen!> Whatever this next thing means, it
    reads

    <\eqnarray*>
      <tformat|<table|<row|<cell|<big|sum><rsub|n=1><rsup|N<rsub|j>><enorm|P<rsub|n><rsup|j>w||2>>|<cell|=>|<cell|<big|sum><rsub|n=1><rsup|N<rsub|j>><enorm|P<rsub|n><rsup|j>P<rsub|j>w||2>\<leqslant\><big|sum><rsub|n=1><rsup|N<rsub|j>><enorm|P<rsub|j>w||>>>>>
    </eqnarray*>

    \;
  </proof>
</body>

<\initial>
  <\collection>
    <associate|language|american>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|7|?>>
    <associate|auto-2|<tuple|7.1|?>>
    <associate|eq:l2-approx-quality|<tuple|7.5|?>>
    <associate|eq:muschwa-diff|<tuple|7.1|?>>
    <associate|eq:muschwa-idsum|<tuple|7.2|?>>
    <associate|eq:stable-decomp|<tuple|7.3|?>>
    <associate|eq:strengthened-cs|<tuple|7.4|?>>
    <associate|eq:strong-cs|<tuple|7.4|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|7<space|2spc>Multiplicative
      Schwarz Methods> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|7.1<space|2spc>V-cycle estimate without
      regularity <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>
    </associate>
  </collection>
</auxiliary>