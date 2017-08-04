import { select } from 'd3';

// toggle visibility of dots or ines in the chart for a group
function toggleShape(entry, shape, chart, isSelected) {
  // toggle the dots on the chart
  const opacity = isSelected ? 0 : 1;
  if (entry === 'Select%20All') {
    chart.selectAll(`.${shape}`)
      .style('opacity', opacity);
  } else {
    chart.selectAll(`[id="id${shape}${entry}"]`)
      .style('opacity', opacity);
  }
}

// used to toggle the color of the key in the legend
function toggleColor(entry, shape, c, color, entries) {
  const clickedLegend = select(`[id="id${shape}${entry}"]`);
  const isSelected = (clickedLegend.style('fill') !== 'white');
  const newColor = isSelected ? 'white' : c;
  clickedLegend.style('fill', newColor);
  if (entry === 'Select%20All') {
    for (const e of entries) {
      select(`[id="id${shape}${e}"]`).style('fill', isSelected ? 'white' : color(e));
    }
  }
  return isSelected;
}

export default function appendLegendKey(legend, i, entry, ly, c, color, entries, chart) {
  // line toggle in the legend
  legend.append('rect')
      .attr('id', `idrect${entry}`)
      .attr('class', 'legend')
      .attr('x', 0)
      .attr('y', ly - 2.5)
      .attr('width', 15)
      .attr('height', 5)
      .style('stroke', 'darkGrey')
      .style('fill', 'white')
      .on('click', () => {
        const b = toggleColor(entry, 'rect', c, color, entries);
        toggleShape(entry, 'line', chart, b);
      });
  // dot toggle in the legend
  legend.append('circle')
      .attr('id', `iddot${entry}`)
      .attr('class', 'legend')
      .attr('cx', 30)
      .attr('cy', ly)
      .attr('r', 5)
      .style('stroke', 'darkGrey')
      .style('fill', c)
      .on('click', () => {
        const b = toggleColor(entry, 'dot', c, color, entries);
        toggleShape(entry, 'circle', chart, b);
      });
  // text for key in the legend
  legend.append('text')
      .attr('class', 'legend')
      .attr('x', 40)
      .attr('y', ly + 5)
      .style('font', '10px sans-serif')
      .text(decodeURI(entry));
}

