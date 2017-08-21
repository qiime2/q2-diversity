/* global categories */
/* global metrics */
import { select } from 'd3';
import state from './state';
import { addMetricPicker, addCategoryPicker } from './toolbar';

export default function init() {
  const metric = metrics[0];
  const category = categories[0];
  const body = select('#main');
  // Order matters here.
  body.append('h1').text('Alpha Rarefaction');
  const controlsRow = body.append('div').attr('class', 'controls row');
  const vizRow = body.append('div').attr('class', 'viz row');

  const plotCol = vizRow.append('div')
                    .attr('class', 'col-lg-9')
                    .style('width', '1120px');
  const legendCol = vizRow.append('div')
                      .attr('class', 'col-lg-2')
                      .style('height', '1000px');
  const plotSvg = plotCol.append('svg');
  const chart = plotSvg.append('g').attr('id', 'chart');
  chart.append('g').attr('class', 'x axis');
  chart.append('g').attr('class', 'y axis');
  chart.append('text').attr('class', 'x label');
  chart.append('text').attr('class', 'y label');
  const subChart = plotSvg.append('g').attr('id', 'subChart');
  subChart.append('g').attr('class', 'x axis');
  subChart.append('g').attr('class', 'y axis');
  subChart.append('text').attr('class', 'x label');
  subChart.append('text').attr('class', 'y label');
  const legendTitle = legendCol.append('div')
                        .style('height', '25px')
                        .style('width', '300px')
                        .style('overflow-y', 'hidden')
                        .style('overflow-x', 'auto')
                        .append('svg')
                          .attr('width', '200')
                          .attr('height', '20')
                          .attr('class', 'legendTitleSvg')
                          .append('g')
                            .attr('class', 'legendTitle');
  const legendBox = legendCol.append('div')
                    .style('height', '900px')
                    .style('width', '300px')
                    .style('overflow-y', 'scroll')
                    .style('overflow-x', 'auto')
                    .attr('class', 'legendBoxDiv')
                    .append('svg')
                      .attr('class', 'legendBoxSvg')
                      .append('g');
  state.initialize(metric, category, controlsRow, plotSvg, legendBox, legendTitle);
  addMetricPicker(controlsRow, metrics, metric);
  if (categories.length > 0) {
    addCategoryPicker(controlsRow, categories, category);
  }
}
