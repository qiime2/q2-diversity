import {
  scaleLinear,
  axisBottom,
  axisLeft,
  scaleOrdinal,
  schemeCategory20,
  select,
  line,
} from 'd3';

import { setupXLabel, setupYLabel } from './axis';

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
  console.log(clickedLegend, clickedLegend.style('fill'), newColor, entry);
  if (entry === 'Select%20All') {
    for (const e of entries) {
      select(`[id="id${shape}${e}"]`).style('fill', isSelected ? 'white' : color(e));
    }
  }
  return isSelected;
}

// add a key to the legend
function appendLegendKey(legend, i, entry, ly, c, color, entries, chart) {
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

// re-render chart and legend whenever selection changes
function renderPlot(svg, data, x, y, category, legend, legendTitle) {
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
  const arrGroups = Array.from(setGroups);

  const legendBox = select(legend.node().parentNode);

  legend.selectAll('.legend').remove();
  legendTitle.selectAll('.legend').remove();
  chart.selectAll('.circle').remove();
  chart.selectAll('.line').remove();

  legend.attr('height', arrGroups.length * 20);
  let ly = 0;

  const valueline = line()
    .x(d => x(d[depthIndex]))
    .y(d => y(d[medianIndex]));

  appendLegendKey(legendTitle, 0, 'Select%20All', 10, 'black',
                  color, arrGroups, chart);
  for (const [i, entry] of arrGroups.entries()) {
    ly = (i + 0.5) * 20;
    const subset = points.filter(d => d[groupIndex] === entry)
                    .sort((a, b) => a[depthIndex] - b[depthIndex]);
    subset.lineOpacity = 0;
    subset.dotOpacity = 1;
    appendLegendKey(legend, i + 1, entry, ly, color(entry), color,
                    arrGroups, chart);
    const curColor = color(subset);
    chart.append('path')
        .attr('d', valueline(subset))
        .style('stroke', curColor)
        .style('opacity', subset.lineOpacity)
        .style('fill', 'none')
        .attr('class', 'line')
        .attr('id', `idline${entry}`);
    chart.selectAll('dot')
        .data(subset)
      .enter()
        .append('circle')
          .attr('cx', d => x(d[depthIndex]))
          .attr('cy', d => y(d[medianIndex]))
          .attr('r', 4)
          .style('stroke', curColor)
          .style('opacity', subset.dotOpacity)
          .style('fill', curColor)
          .attr('class', 'circle')
          .attr('id', d => `idcircle${d[groupIndex]}`);
  }
  legendBox.attr('viewBox', `0 0 200 ${ly + 10}`);
}

// re-render chart edges, exis, formatting, etc. when selection changes
export default function render(svg, data, category, legend, legendTitle) {
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

  renderPlot(svg, data, x, y, category, legend, legendTitle);

  svg.attr('width', width + margin.left + margin.right)
    .attr('height', height + margin.bottom + margin.top);
}
