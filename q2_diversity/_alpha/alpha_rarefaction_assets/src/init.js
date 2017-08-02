/* global categories */
/* global metrics */

import { select } from 'd3';

import { state } from './data';

import { addMetricPicker, addCategoryPicker } from './toolbar';


export default function init() {
  const metric = metrics[0];
  const category = categories[0];
  // DOM
  const body = select('#main');
  const plotRow = body.insert('div', ':first-child').attr('class', 'viz row');
  const plotDiv = plotRow.append('div').attr('class', 'col-md-8');
  const controlsRow = plotDiv.append('div').attr('class', 'controls row');
  const svgRow = plotDiv.append('div').attr('class', 'plot row');
  const svgCol = svgRow.append('div').attr('class', 'col-md-8');
  const svg = svgCol.append('svg');
  const chart = svg.append('g');
  body.insert('h1', ':first-child').text('Alpha Rarefaction');
  chart.append('g').attr('class', 'x axis');
  chart.append('g').attr('class', 'y axis');
  chart.append('text').attr('class', 'x label');
  chart.append('text').attr('class', 'y label');
  // legend
  const legend = plotRow.append('div')
    .attr('class', 'col-md-4')
    .attr('id', 'legendContainer')
    .style('height', chart.height)
    .style('border', '2px solid cyan')
    .style('x', '1200px')
    .append('svg')
      .append('g');

  // D3
  state.initialize(metric, category, controlsRow, svg, legend);
  addMetricPicker(controlsRow, metrics, metric);
  if (categories.length > 0) {
    addCategoryPicker(controlsRow, categories, category);
  }

  // STATS
  // stats(body, data);

  // WARNINGS
  // warnings(body, data);
}
