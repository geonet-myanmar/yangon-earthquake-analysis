/* ── Navigation: active link on scroll ───────────────────────────────────── */
(function () {
  const navbar = document.getElementById('navbar');
  const links  = Array.from(document.querySelectorAll('.nav-links a[href^="#"]'));
  const sections = links.map(a => document.querySelector(a.getAttribute('href'))).filter(Boolean);

  function onScroll() {
    const scrollY = window.scrollY + 80;
    let active = sections[0];
    sections.forEach(s => { if (s.offsetTop <= scrollY) active = s; });
    links.forEach(a => {
      a.classList.toggle('active', a.getAttribute('href') === '#' + active.id);
    });
  }
  window.addEventListener('scroll', onScroll, { passive: true });
  onScroll();
})();

/* ── Lightbox ────────────────────────────────────────────────────────────── */
(function () {
  const lb      = document.getElementById('lightbox');
  const lbImg   = document.getElementById('lightbox-img');
  const lbCap   = document.getElementById('lightbox-caption');

  document.querySelectorAll('.figure-full').forEach(fig => {
    fig.addEventListener('click', () => {
      const img = fig.querySelector('img');
      const cap = fig.querySelector('.figure-caption');
      lbImg.src = img.src;
      lbCap.textContent = cap ? cap.textContent.trim() : '';
      lb.classList.add('open');
      document.body.style.overflow = 'hidden';
    });
  });

  function close() {
    lb.classList.remove('open');
    document.body.style.overflow = '';
  }
  lb.addEventListener('click', close);
  document.addEventListener('keydown', e => { if (e.key === 'Escape') close(); });
})();

/* ── Sortable table ──────────────────────────────────────────────────────── */
(function () {
  const table = document.getElementById('stations-table');
  if (!table) return;

  let sortCol = -1, sortDir = 1;

  table.querySelectorAll('thead th').forEach((th, colIdx) => {
    th.addEventListener('click', () => {
      const tbody = table.querySelector('tbody');
      const rows  = Array.from(tbody.querySelectorAll('tr'));

      if (sortCol === colIdx) {
        sortDir *= -1;
      } else {
        sortDir = 1;
        sortCol = colIdx;
      }

      table.querySelectorAll('thead th').forEach(t => t.classList.remove('sort-asc', 'sort-desc'));
      th.classList.add(sortDir === 1 ? 'sort-asc' : 'sort-desc');

      rows.sort((a, b) => {
        const av = a.cells[colIdx]?.dataset.val ?? a.cells[colIdx]?.textContent.trim() ?? '';
        const bv = b.cells[colIdx]?.dataset.val ?? b.cells[colIdx]?.textContent.trim() ?? '';
        const an = parseFloat(av), bn = parseFloat(bv);
        if (!isNaN(an) && !isNaN(bn)) return (an - bn) * sortDir;
        return av.localeCompare(bv) * sortDir;
      });

      rows.forEach(r => tbody.appendChild(r));
    });
  });
})();

/* ── Animate stat counters on scroll ─────────────────────────────────────── */
(function () {
  const counters = document.querySelectorAll('[data-count]');
  if (!counters.length) return;

  const io = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (!entry.isIntersecting) return;
      const el    = entry.target;
      const end   = parseFloat(el.dataset.count);
      const isInt = Number.isInteger(end);
      const dur   = 900;
      const start = performance.now();

      function step(now) {
        const t   = Math.min((now - start) / dur, 1);
        const val = end * (t < 0.5 ? 2*t*t : -1+(4-2*t)*t);
        el.textContent = isInt ? Math.round(val).toLocaleString()
                                : val.toFixed(2);
        if (t < 1) requestAnimationFrame(step);
        else el.textContent = isInt ? end.toLocaleString() : end.toFixed(2);
      }
      requestAnimationFrame(step);
      io.unobserve(el);
    });
  }, { threshold: 0.5 });

  counters.forEach(c => io.observe(c));
})();
