# =============================================================================
# PLOT POPUP VIEWER MODULE
# =============================================================================
# Provides fullscreen plot viewing with drawing tools and export for any plotly plot.
# Works via client-side JS: click expand icon on any plotly plot to open a
# modal that re-renders the plot with drawing tools and export controls.

# -----------------------------------------------------------------------------
# JS + CSS for automatic plotly plot expansion
# -----------------------------------------------------------------------------

plotPopupJS <- function() {
  tagList(
    # CSS moved to www/app.css

    tags$script(HTML("
      (function() {
        'use strict';

        // --- Expand button injection ---
        // Only add to plots that have finished rendering (have .data property).
        function addExpandBtn(plotDiv) {
          try {
            // Skip if already wrapped or if plot has no data yet
            if (!plotDiv || !plotDiv.data || !plotDiv.layout) return;
            var wrapper = plotDiv.closest('.plotly-expand-wrapper');
            if (wrapper) return; // already processed

            var parent = plotDiv.parentElement;
            if (!parent) return;

            // Wrap the plot's parent in a relative container for the button
            var wrap = document.createElement('div');
            wrap.className = 'plotly-expand-wrapper';
            parent.insertBefore(wrap, plotDiv);
            wrap.appendChild(plotDiv);

            var btn = document.createElement('div');
            btn.className = 'plotly-expand-btn';
            btn.innerHTML = '<i class=\"fa fa-expand\"></i> Expand';
            btn.title = 'Open in fullscreen viewer';
            btn.addEventListener('click', function(e) {
              e.stopPropagation();
              e.preventDefault();
              openPopup(plotDiv);
            });
            wrap.appendChild(btn);
          } catch(err) {
            // Silently ignore — don't break the app
          }
        }

        // Debounced scan: avoid hammering the DOM during rapid updates
        var _scanTimer = null;
        function scheduleScan() {
          if (_scanTimer) return;
          _scanTimer = setTimeout(function() {
            _scanTimer = null;
            document.querySelectorAll('.js-plotly-plot').forEach(addExpandBtn);
          }, 500);
        }

        // Observe DOM for new plotly plots (debounced)
        var observer = new MutationObserver(scheduleScan);
        if (document.body) {
          observer.observe(document.body, { childList: true, subtree: true });
        } else {
          document.addEventListener('DOMContentLoaded', function() {
            observer.observe(document.body, { childList: true, subtree: true });
          });
        }

        // Periodic fallback scan (every 5s, lightweight)
        setInterval(function() {
          document.querySelectorAll('.js-plotly-plot').forEach(addExpandBtn);
        }, 5000);

        // Close on Escape
        document.addEventListener('keydown', function(e) {
          if (e.key === 'Escape') closePopup();
        });

        // --- Modal management ---
        var _currentData = null;
        var _currentLayout = null;

        function openPopup(plotDiv) {
          // Extract plotly data from the DOM element
          try {
            if (!plotDiv.data || !plotDiv.layout) {
              console.warn('Plot not ready for expansion');
              return;
            }
            _currentData = JSON.parse(JSON.stringify(plotDiv.data));
            _currentLayout = JSON.parse(JSON.stringify(plotDiv.layout));
          } catch(err) {
            console.warn('Could not extract plot data:', err);
            return;
          }

          var overlay = document.getElementById('plotPopupOverlay');
          if (!overlay) {
            createModal();
            overlay = document.getElementById('plotPopupOverlay');
          }

          overlay.style.display = 'block';
          renderPopupPlot();
        }

        function closePopup() {
          var overlay = document.getElementById('plotPopupOverlay');
          if (overlay) overlay.style.display = 'none';
          var container = document.getElementById('plotPopupContainer');
          if (container) {
            try { Plotly.purge(container); } catch(e) {}
          }
        }

        function renderPopupPlot() {
          if (!_currentData || !_currentLayout) return;
          var container = document.getElementById('plotPopupContainer');
          if (!container) return;
          var layout = Object.assign({}, _currentLayout);

          // Adjust margins for fullscreen
          layout.margin = { l: 60, r: 30, t: 50, b: 60 };
          layout.autosize = true;

          // Drawing tools config
          var config = {
            displayModeBar: true,
            modeBarButtonsToAdd: [
              'drawline', 'drawcircle', 'drawrect', 'eraseshape'
            ],
            toImageButtonOptions: {
              format: 'png',
              height: 1200,
              width: 1600,
              scale: 2
            },
            scrollZoom: true,
            displaylogo: false
          };

          Plotly.newPlot(container, _currentData, layout, config).then(function() {
            Plotly.Plots.resize(container);
          }).catch(function(err) {
            console.warn('Popup plot render error:', err);
          });
        }

        // Export function
        window._popupExport = function(format) {
          var container = document.getElementById('plotPopupContainer');
          if (!container || !container.data) return;
          var dpiSelect = document.getElementById('popupDPI');
          var scale = dpiSelect ? parseFloat(dpiSelect.value) : 2;
          var widthInput = document.getElementById('popupWidth');
          var heightInput = document.getElementById('popupHeight');
          var w = widthInput ? parseInt(widthInput.value) : 1600;
          var h = heightInput ? parseInt(heightInput.value) : 1000;

          Plotly.toImage(container, {
            format: format,
            width: w,
            height: h,
            scale: scale
          }).then(function(url) {
            var a = document.createElement('a');
            a.href = url;
            a.download = 'plot_export_' + new Date().toISOString().slice(0,10) + '.' + format;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
          });
        };

        // Clear shapes
        window._popupClearShapes = function() {
          var container = document.getElementById('plotPopupContainer');
          if (container) {
            Plotly.relayout(container, { shapes: [] });
          }
        };

        // Resize on window resize
        window.addEventListener('resize', function() {
          var container = document.getElementById('plotPopupContainer');
          var overlay = document.getElementById('plotPopupOverlay');
          if (container && overlay && overlay.style.display === 'block') {
            Plotly.Plots.resize(container);
          }
        });

        // --- Build modal DOM ---
        function createModal() {
          var overlay = document.createElement('div');
          overlay.id = 'plotPopupOverlay';
          overlay.addEventListener('click', function(e) {
            if (e.target === overlay) closePopup();
          });

          var modal = document.createElement('div');
          modal.id = 'plotPopupModal';

          // Toolbar
          var toolbar = document.createElement('div');
          toolbar.id = 'plotPopupToolbar';
          toolbar.innerHTML = [
            '<div class=\"toolbar-group\">',
              '<label>Width:</label>',
              '<input id=\"popupWidth\" type=\"number\" value=\"1600\" min=\"400\" max=\"4000\" step=\"100\" style=\"width:70px;\">',
            '</div>',
            '<div class=\"toolbar-group\">',
              '<label>Height:</label>',
              '<input id=\"popupHeight\" type=\"number\" value=\"1000\" min=\"300\" max=\"3000\" step=\"100\" style=\"width:70px;\">',
            '</div>',
            '<div class=\"toolbar-group\">',
              '<label>Scale:</label>',
              '<select id=\"popupDPI\">',
                '<option value=\"1\">1x</option>',
                '<option value=\"2\" selected>2x (Recommended)</option>',
                '<option value=\"3\">3x</option>',
                '<option value=\"4\">4x (Publication)</option>',
              '</select>',
            '</div>',
            '<div class=\"toolbar-sep\"></div>',
            '<button onclick=\"_popupExport(\\x27png\\x27)\"><i class=\"fa fa-image\"></i> PNG</button>',
            '<button onclick=\"_popupExport(\\x27svg\\x27)\"><i class=\"fa fa-vector-square\"></i> SVG</button>',
            '<div class=\"toolbar-sep\"></div>',
            '<button class=\"btn-danger\" onclick=\"_popupClearShapes()\"><i class=\"fa fa-eraser\"></i> Clear Shapes</button>',
            '<div style=\"flex:1;\"></div>',
            '<button class=\"btn-close-popup\" onclick=\"document.getElementById(\\x27plotPopupOverlay\\x27).style.display=\\x27none\\x27; try{Plotly.purge(document.getElementById(\\x27plotPopupContainer\\x27));}catch(e){}\">',
              '<i class=\"fa fa-times\"></i> Close (Esc)</button>'
          ].join('');

          var container = document.createElement('div');
          container.id = 'plotPopupContainer';

          modal.appendChild(toolbar);
          modal.appendChild(container);
          overlay.appendChild(modal);
          document.body.appendChild(overlay);
        }

      })();
    "))
  )
}
