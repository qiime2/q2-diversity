import {
  scaleLinear,
  axisBottom,
  axisLeft,
  scaleOrdinal,
  schemeCategory20,
  select,
} from 'd3';

import { setupXLabel, setupYLabel } from './axis';

// toggle visibility of dots in the chart for a group
function toggleDots(chart, isSelected, value) {
  // toggle the dots on the chart
  const opacity = isSelected ? 0 : 1;
  console.log('selectAll(', `.dot${value.replace(' ', '-')}`, ') ',
              'results: ', chart.selectAll(`.dot${value.replace(' ', '-')}`),
              'from data: ', chart.selectAll('circle'));
  chart.selectAll(`.dot${value.replace(' ', '-')}`)
    .style('opacity', opacity);
}

// toggle visibility of a line in the chart for a group
function toggleLine(entry, isSelected) {
  // toggle the line in the entry
  console.log(entry, isSelected);
}

// used to toggle the color of the item in the legend
function toggleColor(entry, shape, c, color, entries) {
  const clickedLegend = select(`#id${shape}${entry.replace(' ', '-')}`);
  const isSelected = (clickedLegend.style('fill') !== 'white');
  const newColor = isSelected ? 'white' : c;
  clickedLegend.style('fill', newColor);
  if (entry === 'Select All') {
    for (const e of entries) {
      select(`#id${shape}${e.replace(' ', '-')}`).style('fill', isSelected ? 'white' : color(e));
    }
  }
  return isSelected;
}

// add a key to the legend
function appendLegendKey(legend, i, entry, ly, c, color, entries, chart) {
  // line toggle in the legend
  legend.append('rect')
      .attr('id', `idrect${entry.replace(' ', '-')}`)
      .attr('class', 'legend')
      .attr('x', 0)
      .attr('y', ly - 2.5)
      .attr('width', 15)
      .attr('height', 5)
      .style('stroke', 'darkGrey')
      .style('fill', 'white')
      .on('click', () => {
        const b = toggleColor(entry, 'rect', c, color, entries);
        toggleLine(entry, b);
      });
  // dot toggle in the legend
  legend.append('circle')
      .attr('id', `idcircle${entry.replace(' ', '-')}`)
      .attr('class', 'legend')
      .attr('cx', 30)
      .attr('cy', ly)
      .attr('r', 5)
      .style('stroke', 'darkGrey')
      .style('fill', c)
      .on('click', () => {
        const b = toggleColor(entry, 'circle', c, color, entries);
        toggleDots(chart, b, entry);
      });
  // text for key in the legend
  legend.append('text')
      .attr('class', 'legend')
      .attr('x', 40)
      .attr('y', ly + 5)
      .style('font', '10px sans-serif')
      .text(entry);
}

// re-render chart and legend whenever selection changes
function renderPlot(svg, data, x, y, category, legend) {
  const chart = svg.select('g');

  const depthIndex = data.data.columns.indexOf('depth');
  const medianIndex = data.data.columns.indexOf('median');
  let groupIndex = data.data.columns.indexOf('sample-id');
  if (groupIndex === -1) {
    groupIndex = data.data.columns.indexOf(category);
  }
  const points = [data.data.data.sort((a, b) => a[depthIndex] - b[depthIndex])][0];
  const setGroups = new Set(Array.from(points, d => d[groupIndex]));
  const color = scaleOrdinal(schemeCategory20)
    .domain(setGroups);

  chart.selectAll('circle').remove();
  chart.selectAll('dot')
      .data(points)
    .enter()
      .append('circle')
        .attr('cx', d => x(d[depthIndex]))
        .attr('cy', d => y(d[medianIndex]))
        .attr('r', 4)
        .style('stroke', d => color(d[groupIndex]))
        .style('fill', d => color(d[groupIndex]))
        .attr('class', d => `dot${d[groupIndex].replace(' ', '-')}`);

  legend.selectAll('.legend').remove();
  const arrGroups = Array.from(setGroups);
  legend.attr('height', arrGroups.length * 20);
  let ly = 0;
  const legendBox = select(legend.node().parentNode);
  appendLegendKey(legendBox, 0, 'Select All', 10, 'black',
                  color, arrGroups, chart);
  for (const [i, entry] of arrGroups.entries()) {
    ly = (i + 1.5) * 20;
    const c = color(entry);
    appendLegendKey(legend, i + 1, entry, ly, c, color,
                    arrGroups, chart);
  }
  legendBox.attr('viewBox', `0 0 200 ${ly + 10}`);
}

// re-render chart edges, exis, formatting, etc. when selection changes
export default function render(svg, data, category, legend) {
  const height = 400;
  const width = 1000;
  const margin = { top: 20, left: 70, right: 50, bottom: 50 };
  const chart = svg.select('g');

  const { xAxisLabel, yAxisLabel, minX, maxX, minY, maxY } = data;

  const xAxis = axisBottom();
  const yAxis = axisLeft();

  let pad = (maxX - minX) * 0.03;
  if (Number.isInteger(minX) && Number.isInteger(maxX)) {
    pad = Math.max(Math.round(pad), 1);
    const between = Math.max(3, (maxX - minX) + (2 * pad));
    xAxis.ticks(Math.min(between, 12), 'd');
  }

  const x = scaleLinear().domain([minX - pad, maxX + pad]).range([0, width]).nice();
  const y = scaleLinear().domain([minY, maxY]).range([height, 0]).nice();

  xAxis.scale(x);
  yAxis.scale(y);

  chart.attr('transform', `translate(${margin.left},${margin.top})`);

  setupXLabel(svg, width, height, xAxisLabel, xAxis);
  setupYLabel(svg, height, yAxisLabel, yAxis);

  renderPlot(svg, data, x, y, category, legend);

  svg.attr('width', width + margin.left + margin.right)
    .attr('height', height + margin.bottom + margin.top);
}
